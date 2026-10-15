#include "GaussSeidelNoCopySmoother.H"
#include "PrecisionAdaptor.H"

namespace Foam
{
    defineTypeNameAndDebug(GaussSeidelNoCopySmoother, 0);

    lduMatrix::smoother::addsymMatrixConstructorToTable<GaussSeidelNoCopySmoother>
        addGaussSeidelNoCopySmootherSymMatrixConstructorToTable_;

    lduMatrix::smoother::addasymMatrixConstructorToTable<GaussSeidelNoCopySmoother>
        addGaussSeidelNoCopySmootherAsymMatrixConstructorToTable_;
}


Foam::GaussSeidelNoCopySmoother::GaussSeidelNoCopySmoother
(
    const word& fieldName,
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const FieldField<Field, scalar>& interfaceIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const dictionary& solverControls
)
:
    lduMatrix::smoother
    (
        fieldName,
        matrix,
        interfaceBouCoeffs,
        interfaceIntCoeffs,
        interfaces
    )
{
    const label nCells = matrix.diag().size();
    const labelUList& u = matrix.lduAddr().upperAddr();

    boolList isBoundary(nCells, false);
    forAll(interfaces, patchi)
    {
        if (interfaces.set(patchi))
        {
            for (const label celli : interfaces[patchi].interface().faceCells())
            {
                isBoundary[celli] = true;
            }
        }
    }

    // Faces are visited in index order by the sweep (rows in order,
    // ownerStartAddr is monotonic), so the first face with a given upper
    // cell makes that cell's first update
    boolList written(nCells, false);
    firstWrite_.resize_nocopy(u.size());
    forAll(u, facei)
    {
        const label celli = u[facei];
        firstWrite_[facei] = !written[celli] && !isBoundary[celli];
        written[celli] = true;
    }

    readSource_.resize_nocopy(nCells);
    DynamicList<label> bCells;
    for (label celli = 0; celli < nCells; ++celli)
    {
        readSource_[celli] = !written[celli] && !isBoundary[celli];
        if (isBoundary[celli]) bCells.push_back(celli);
    }
    boundaryCells_.transfer(bCells);
}


void Foam::GaussSeidelNoCopySmoother::sweep
(
    solveScalarField& psi,
    const solveScalarField& source,
    const direction cmpt,
    const label nSweeps
) const
{
    solveScalar* __restrict__ psiPtr = psi.begin();
    const solveScalar* const __restrict__ srcPtr = source.begin();

    const label nCells = psi.size();

    solveScalarField& bPrime = matrix_.work(nCells);
    solveScalar* __restrict__ bPrimePtr = bPrime.begin();

    const scalar* const __restrict__ diagPtr = matrix_.diag().begin();
    const scalar* const __restrict__ upperPtr = matrix_.upper().begin();
    const scalar* const __restrict__ lowerPtr = matrix_.lower().begin();

    const label* const __restrict__ uPtr =
        matrix_.lduAddr().upperAddr().begin();
    const label* const __restrict__ ownStartPtr =
        matrix_.lduAddr().ownerStartAddr().begin();

    const bool* const __restrict__ firstPtr = firstWrite_.begin();
    const bool* const __restrict__ readSrcPtr = readSource_.begin();

    for (label sweep=0; sweep<nSweeps; sweep++)
    {
        // Only the cells the coupled-interface update writes into
        for (const label celli : boundaryCells_)
        {
            bPrimePtr[celli] = srcPtr[celli];
        }

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

        solveScalar psii;
        label fStart;
        label fEnd = ownStartPtr[0];

        for (label celli=0; celli<nCells; celli++)
        {
            fStart = fEnd;
            fEnd = ownStartPtr[celli + 1];

            psii = readSrcPtr[celli] ? srcPtr[celli] : bPrimePtr[celli];

            for (label facei=fStart; facei<fEnd; facei++)
            {
                psii -= upperPtr[facei]*psiPtr[uPtr[facei]];
            }

            psii /= diagPtr[celli];

            for (label facei=fStart; facei<fEnd; facei++)
            {
                const label nbr = uPtr[facei];
                if (firstPtr[facei])
                {
                    bPrimePtr[nbr] = srcPtr[nbr] - lowerPtr[facei]*psii;
                }
                else
                {
                    bPrimePtr[nbr] -= lowerPtr[facei]*psii;
                }
            }

            psiPtr[celli] = psii;
        }
    }
}


void Foam::GaussSeidelNoCopySmoother::smooth
(
    solveScalarField& psi,
    const scalarField& source,
    const direction cmpt,
    const label nSweeps
) const
{
    sweep(psi, ConstPrecisionAdaptor<solveScalar, scalar>(source), cmpt, nSweeps);
}


void Foam::GaussSeidelNoCopySmoother::scalarSmooth
(
    solveScalarField& psi,
    const solveScalarField& source,
    const direction cmpt,
    const label nSweeps
) const
{
    sweep(psi, source, cmpt, nSweeps);
}
