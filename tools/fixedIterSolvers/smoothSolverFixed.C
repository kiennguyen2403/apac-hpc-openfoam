#include "smoothSolverFixed.H"
#include "fixedIterControl.H"
#include "PrecisionAdaptor.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    defineTypeNameAndDebug(smoothSolverFixed, 0);

    lduMatrix::solver::addsymMatrixConstructorToTable<smoothSolverFixed>
        addsmoothSolverFixedSymMatrixConstructorToTable_;

    lduMatrix::solver::addasymMatrixConstructorToTable<smoothSolverFixed>
        addsmoothSolverFixedAsymMatrixConstructorToTable_;
}


Foam::smoothSolverFixed::smoothSolverFixed
(
    const word& fieldName,
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const FieldField<Field, scalar>& interfaceIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const dictionary& solverControls
)
:
    smoothSolver
    (
        fieldName,
        matrix,
        interfaceBouCoeffs,
        interfaceIntCoeffs,
        interfaces,
        solverControls
    )
{}


Foam::solverPerformance Foam::smoothSolverFixed::solve
(
    scalarField& psi_s,
    const scalarField& source,
    const direction cmpt
) const
{
    if (nSweeps_ < 0 || !convergenceImpossible(tolerance_, relTol_, log_))
    {
        return smoothSolver::solve(psi_s, source, cmpt);
    }

    // As smoothSolver::solve with nSweeps_ >= 0, minus the residual norm
    // after every sweep but the last (checkConvergence is always false).
    PrecisionAdaptor<solveScalar, scalar> tpsi(psi_s);
    solveScalarField& psi = tpsi.ref();

    solverPerformance solverPerf(smoothSolver::typeName, fieldName_);

    solveScalar normFactor = 0;
    solveScalarField residual;

    ConstPrecisionAdaptor<solveScalar, scalar> tsource(source);

    {
        solveScalarField Apsi(psi.size());
        solveScalarField temp(psi.size());

        matrix_.Amul(Apsi, psi, interfaceBouCoeffs_, interfaces_, cmpt);

        normFactor = this->normFactor(psi, tsource(), Apsi, temp);

        residual = tsource() - Apsi;

        matrix().setResidualField
        (
            ConstPrecisionAdaptor<scalar, solveScalar>(residual)(),
            fieldName_,
            true
        );

        solverPerf.initialResidual() =
            gSumMag(residual, matrix().mesh().comm())/normFactor;
        solverPerf.finalResidual() = solverPerf.initialResidual();
    }

    {
        addProfiling(solve, "lduMatrix::smoother.", fieldName_);

        autoPtr<lduMatrix::smoother> smootherPtr = lduMatrix::smoother::New
        (
            fieldName_,
            matrix_,
            interfaceBouCoeffs_,
            interfaceIntCoeffs_,
            interfaces_,
            controlDict_
        );

        if (nSweeps_ > 0 && smootherPtr->type() == "GaussSeidel")
        {
            // v2: GaussSeidel does n sweeps per call exactly as n calls of
            // one sweep (per sweep: bPrime = source, interfaces, sweep), so
            // all the sweeps of the loop below go into one call. Bit-identical,
            // and lets a smoother fold the per-sweep copy (tools/stepReuse).
            label nCalls = 0;
            do
            {
                ++nCalls;
            } while
            (
                (solverPerf.nIterations() += nSweeps_) < maxIter_
             || solverPerf.nIterations() < minIter_
            );

            smootherPtr->smooth(psi, source, cmpt, nCalls*nSweeps_);
        }
        else
        {
            do
            {
                smootherPtr->smooth(psi, source, cmpt, nSweeps_);
            } while
            (
                (solverPerf.nIterations() += nSweeps_) < maxIter_
             || solverPerf.nIterations() < minIter_
            );
        }

        residual =
            matrix_.residual(psi, source, interfaceBouCoeffs_, interfaces_, cmpt);

        solverPerf.finalResidual() =
            gSumMag(residual, matrix().mesh().comm())/normFactor;
    }

    matrix().setResidualField
    (
        ConstPrecisionAdaptor<scalar, solveScalar>(residual)(),
        fieldName_,
        false
    );

    return solverPerf;
}
