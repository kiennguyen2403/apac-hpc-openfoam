/*---------------------------------------------------------------------------*\
    compactAddr: LD_PRELOAD replacements of the two linear-solver kernels that
    stream the matrix addressing on every sweep:

      GaussSeidelSmoother::smooth (static worker; also keeps libstepReuse's
                                   folded per-sweep source copy)
      lduMatrix::Amul

    The solver time loop is memory-bandwidth bound (job 182171696: 72% of the
    STREAM triad bandwidth, 72% of cycles with >= 4 reads outstanding), so the
    addressing is stored in a compact form, built once per mesh level and
    cached:
      neighbour offset  u - l  as uint16 (instead of the int32 upperAddr)
      faces per owner   as uint8  (instead of the int32 ownerStartAddr)
    Amul no longer reads lowerAddr either: faces are ordered by owner, so the
    owner of a face is the cell whose face range is being walked.

    Gauss-Seidel: per cell 1 B instead of 4 B, per face 2 B instead of 4 B.
    Amul:         per cell 1 B extra, per face 2 B instead of 8 B.

    Same floating-point operations in the same order (faces are visited in
    increasing face index, exactly as the stock loops do): bit-identical.

    A level whose offsets exceed 65535, or whose faces are not owner-ordered,
    runs the next definition (libstepReuse's folded Gauss-Seidel / stock
    Amul). STEP_REUSE_OFF=compact switches both replacements off;
    STEP_REUSE_REPORT=1 prints one line per path on the master.

    LD_PRELOAD order: this library before libstepReuse.
\*---------------------------------------------------------------------------*/

#include "GaussSeidelSmoother.H"
#include "lduMatrix.H"
#include "lduAddressing.H"
#include "UPstream.H"
#include "compactAddressing.H"

#include <dlfcn.h>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <set>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <vector>

namespace
{

void report(const std::string& msg)
{
    static std::set<std::string> done;
    if (!std::getenv("STEP_REUSE_REPORT") || !Foam::UPstream::master())
    {
        return;
    }
    if (done.insert(msg).second)
    {
        std::cerr << "compactAddr: " << msg << std::endl;
    }
}

void* nextSym(const char* name)
{
    void* p = ::dlsym(RTLD_NEXT, name);
    if (!p)
    {
        std::fprintf(stderr, "compactAddr: no next definition of %s\n", name);
        std::abort();
    }
    return p;
}

} // End anonymous namespace

using Foam::compactAddr::Compact;
using Foam::compactAddr::compactOf;
using Foam::compactAddr::off;


void Foam::GaussSeidelSmoother::smooth
(
    const word& fieldName_,
    solveScalarField& psi,
    const lduMatrix& matrix_,
    const solveScalarField& source,
    const FieldField<Field, scalar>& interfaceBouCoeffs_,
    const lduInterfaceFieldPtrsList& interfaces_,
    const direction cmpt,
    const label nSweeps
)
{
    const Compact* cp = off() ? nullptr : compactOf(matrix_.lduAddr());

    if (!cp)
    {
        typedef void (*smoothFn)
        (
            const word&, solveScalarField&, const lduMatrix&,
            const solveScalarField&, const FieldField<Field, scalar>&,
            const lduInterfaceFieldPtrsList&, const direction, const label
        );
        static const smoothFn next = reinterpret_cast<smoothFn>
        (
            nextSym
            (
                "_ZN4Foam19GaussSeidelSmoother6smoothERKNS_4wordERNS_5FieldIdEE"
                "RKNS_9lduMatrixERKS5_RKNS_10FieldFieldIS4_dEE"
                "RKNS_8UPtrListIKNS_17lduInterfaceFieldEEEhi"
            )
        );
        if (!off())
        {
            report("GaussSeidel: level does not fit 16-bit offsets, next definition used");
        }
        next
        (
            fieldName_, psi, matrix_, source, interfaceBouCoeffs_,
            interfaces_, cmpt, nSweeps
        );
        return;
    }

    report("GaussSeidel: 16-bit neighbour offsets, folded source copy");

    solveScalar* __restrict__ psiPtr = psi.begin();

    const label nCells = psi.size();

    solveScalarField& bPrime = matrix_.work(nCells);

    solveScalar* __restrict__ bPrimePtr = bPrime.begin();
    const solveScalar* const __restrict__ sourcePtr = source.begin();

    const scalar* const __restrict__ diagPtr = matrix_.diag().begin();
    const scalar* const __restrict__ upperPtr = matrix_.upper().begin();
    const scalar* const __restrict__ lowerPtr = matrix_.lower().begin();

    const std::uint16_t* const __restrict__ offPtr = cp->off.data();
    const std::uint8_t* const __restrict__ cntPtr = cp->cnt.data();

    // Sweep over the cells (stock order); with the fold,
    // bPrime[celli] = source[celli] once its value for this sweep is read
    auto sweepCells = [&](auto foldTag)
    {
        solveScalar psii;
        label fStart;
        label fEnd = 0;

        for (label celli=0; celli<nCells; celli++)
        {
            // Start and end of this row
            fStart = fEnd;
            fEnd = fStart + cntPtr[celli];

            // Get the accumulated neighbour side
            psii = bPrimePtr[celli];

            if constexpr (decltype(foldTag)::value)
            {
                bPrimePtr[celli] = sourcePtr[celli];
            }

            // Accumulate the owner product side
            for (label facei=fStart; facei<fEnd; facei++)
            {
                psii -= upperPtr[facei]*psiPtr[celli + offPtr[facei]];
            }

            // Finish psi for this cell
            psii /= diagPtr[celli];

            // Distribute the neighbour side using psi for this cell
            for (label facei=fStart; facei<fEnd; facei++)
            {
                bPrimePtr[celli + offPtr[facei]] -= lowerPtr[facei]*psii;
            }

            psiPtr[celli] = psii;
        }
    };

    bPrime = source;

    for (label sweep=0; sweep<nSweeps; sweep++)
    {
        const label startRequest = UPstream::nRequests();

        matrix_.initMatrixInterfaces
        (
            false,
            interfaceBouCoeffs_,
            interfaces_,
            psi,
            bPrime,
            cmpt
        );

        matrix_.updateMatrixInterfaces
        (
            false,
            interfaceBouCoeffs_,
            interfaces_,
            psi,
            bPrime,
            cmpt,
            startRequest
        );

        if (sweep + 1 < nSweeps)
        {
            sweepCells(std::true_type());
        }
        else
        {
            sweepCells(std::false_type());
        }
    }
}


void Foam::lduMatrix::Amul
(
    solveScalarField& Apsi,
    const tmp<solveScalarField>& tpsi,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const direction cmpt
) const
{
    const Compact* cp =
        (off() || hasLowerCSR()) ? nullptr : compactOf(lduAddr());

    if (!cp)
    {
        typedef void (*amulFn)
        (
            const lduMatrix*, solveScalarField&, const tmp<solveScalarField>&,
            const FieldField<Field, scalar>&, const lduInterfaceFieldPtrsList&,
            const direction
        );
        static const amulFn next = reinterpret_cast<amulFn>
        (
            nextSym
            (
                "_ZNK4Foam9lduMatrix4AmulERNS_5FieldIdEERKNS_3tmpIS2_EERKNS_"
                "10FieldFieldIS1_dEERKNS_8UPtrListIKNS_17lduInterfaceFieldEEEh"
            )
        );
        if (!off())
        {
            report("Amul: level does not fit 16-bit offsets, stock used");
        }
        next(this, Apsi, tpsi, interfaceBouCoeffs, interfaces, cmpt);
        return;
    }

    report("Amul: 16-bit neighbour offsets, owner from the face ranges");

    solveScalar* __restrict__ ApsiPtr = Apsi.begin();

    const solveScalarField& psi = tpsi();
    const solveScalar* const __restrict__ psiPtr = psi.begin();

    const scalar* const __restrict__ diagPtr = diag().begin();
    const scalar* const __restrict__ upperPtr = upper().begin();
    const scalar* const __restrict__ lowerPtr = lower().begin();

    const std::uint16_t* const __restrict__ offPtr = cp->off.data();
    const std::uint8_t* const __restrict__ cntPtr = cp->cnt.data();

    const label startRequest = UPstream::nRequests();

    // Initialise the update of interfaced interfaces
    initMatrixInterfaces
    (
        true,
        interfaceBouCoeffs,
        interfaces,
        psi,
        Apsi,
        cmpt
    );

    const label nCells = diag().size();

    for (label cell=0; cell<nCells; cell++)
    {
        ApsiPtr[cell] = diagPtr[cell]*psiPtr[cell];
    }

    // Faces in increasing index: owner = cell, neighbour = cell + offset
    label face = 0;
    for (label cell=0; cell<nCells; cell++)
    {
        const label fEnd = face + cntPtr[cell];
        for (; face<fEnd; face++)
        {
            const label nei = cell + offPtr[face];
            ApsiPtr[nei] += lowerPtr[face]*psiPtr[cell];
            ApsiPtr[cell] += upperPtr[face]*psiPtr[nei];
        }
    }

    // Update interface interfaces
    updateMatrixInterfaces
    (
        true,
        interfaceBouCoeffs,
        interfaces,
        psi,
        Apsi,
        cmpt,
        startRequest
    );

    tpsi.clear();
}


// ************************************************************************* //
