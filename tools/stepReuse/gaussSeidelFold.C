/*---------------------------------------------------------------------------*\
    Replacement for GaussSeidelSmoother::smooth (the static worker that both
    virtual smooth() and scalarSmooth() call through the PLT).

    Stock: every sweep starts with bPrime = source, a separate streaming pass.
    Here that copy, for every sweep after the first, happens inside the
    previous sweep, right after bPrime[celli] has been read for the last time
    in that sweep: later cells only write bPrime[j] for j > their own index,
    so bPrime[celli] is not touched again. The interface update of the next
    sweep then starts from bPrime == source, as in the stock code. Same values
    and the same floating-point operations in the same order: bit-identical.

    One call with a single sweep has nothing to fold and runs the stock code.
\*---------------------------------------------------------------------------*/

#include "GaussSeidelSmoother.H"
#include "stepReuse.H"

#include <type_traits>

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
    if (nSweeps < 2 || stepReuse::off("gs"))
    {
        typedef void (*smoothFn)
        (
            const word&,
            solveScalarField&,
            const lduMatrix&,
            const solveScalarField&,
            const FieldField<Field, scalar>&,
            const lduInterfaceFieldPtrsList&,
            const direction,
            const label
        );

        static const smoothFn stock = reinterpret_cast<smoothFn>
        (
            stepReuse::next
            (
                "_ZN4Foam19GaussSeidelSmoother6smoothERKNS_4wordERNS_5FieldIdEE"
                "RKNS_9lduMatrixERKS5_RKNS_10FieldFieldIS4_dEE"
                "RKNS_8UPtrListIKNS_17lduInterfaceFieldEEEhi"
            )
        );

        stock
        (
            fieldName_, psi, matrix_, source, interfaceBouCoeffs_,
            interfaces_, cmpt, nSweeps
        );
        return;
    }

    stepReuse::report("gs", "GaussSeidel: folded per-sweep source copy");

    solveScalar* __restrict__ psiPtr = psi.begin();

    const label nCells = psi.size();

    solveScalarField& bPrime = matrix_.work(nCells);

    solveScalar* __restrict__ bPrimePtr = bPrime.begin();
    const solveScalar* const __restrict__ sourcePtr = source.begin();

    const scalar* const __restrict__ diagPtr = matrix_.diag().begin();
    const scalar* const __restrict__ upperPtr =
        matrix_.upper().begin();
    const scalar* const __restrict__ lowerPtr =
        matrix_.lower().begin();

    const label* const __restrict__ uPtr =
        matrix_.lduAddr().upperAddr().begin();

    const label* const __restrict__ ownStartPtr =
        matrix_.lduAddr().ownerStartAddr().begin();

    // Sweep over the cells; with the fold, bPrime[celli] = source[celli]
    // once its value for this sweep has been read
    auto sweepCells = [&](auto foldTag)
    {
        solveScalar psii;
        label fStart;
        label fEnd = ownStartPtr[0];

        for (label celli=0; celli<nCells; celli++)
        {
            // Start and end of this row
            fStart = fEnd;
            fEnd = ownStartPtr[celli + 1];

            // Get the accumulated neighbour side
            psii = bPrimePtr[celli];

            if constexpr (decltype(foldTag)::value)
            {
                bPrimePtr[celli] = sourcePtr[celli];
            }

            // Accumulate the owner product side
            for (label facei=fStart; facei<fEnd; facei++)
            {
                psii -= upperPtr[facei]*psiPtr[uPtr[facei]];
            }

            // Finish psi for this cell
            psii /= diagPtr[celli];

            // Distribute the neighbour side using psi for this cell
            for (label facei=fStart; facei<fEnd; facei++)
            {
                bPrimePtr[uPtr[facei]] -= lowerPtr[facei]*psii;
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


// ************************************************************************* //
