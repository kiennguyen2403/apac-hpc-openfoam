/*---------------------------------------------------------------------------*\
    LD_PRELOAD replacement for Foam::smoothSolver::solve, the stock solver
    that the shipped fvSolution selects for U, k and omega.

    With fixed-iteration controls (tolerance 0, relTol 0, as in the shipped
    fvSolution.fixedIter) checkConvergence can never return true and prints
    nothing, so the residual after every sweep (a full matrix pass plus a
    global reduction) only feeds a test that always fails. Here it is
    computed once, after the last sweep: same sweeps, same psi, same initial
    and final residuals and the same log line as the stock solver.

    Any other controls (fixedIterControl.H), nSweeps < 0, or
    STEP_REUSE_OFF=fixed run the stock function (next definition, in
    libOpenFOAM).

    The fixed path is the code of smoothSolverFixed v2 (lib_v2/src), which
    was selected by renaming the solver in fvSolution and loading the library
    with controlDict libs; now the shipped fvSolution and controlDict are used
    unchanged and the library is preloaded.

    GaussSeidel does n sweeps per call exactly as n calls of one sweep (per
    sweep: bPrime = source, interfaces, sweep), so all the sweeps go into one
    call, which lets libcompactAddr / libstepReuse fold the per-sweep copy.
\*---------------------------------------------------------------------------*/

#include "smoothSolver.H"
#include "profiling.H"
#include "PrecisionAdaptor.H"
#include "fixedIterControl.H"
#include "fixedIterSolvers.H"


Foam::solverPerformance Foam::smoothSolver::solve
(
    scalarField& psi_s,
    const scalarField& source,
    const direction cmpt
) const
{
    if
    (
        nSweeps_ < 0
     || fixedIterSolvers::off("fixed")
     || !convergenceImpossible(tolerance_, relTol_, log_)
    )
    {
        // Itanium C++ ABI: the hidden return-value address comes first, then
        // this, as for a free function taking the object as first argument
        typedef solverPerformance (*solveFn)
        (
            const smoothSolver*,
            scalarField&,
            const scalarField&,
            const direction
        );
        static const solveFn stock = reinterpret_cast<solveFn>
        (
            fixedIterSolvers::next
            (
                "_ZNK4Foam12smoothSolver5solveERNS_5FieldIdEERKS2_h"
            )
        );

        if (fixedIterSolvers::reporting())
        {
            fixedIterSolvers::report
            (
                "smoothSolver " + fieldName_ + ": stock solve"
            );
        }

        return stock(this, psi_s, source, cmpt);
    }

    if (fixedIterSolvers::reporting())
    {
        fixedIterSolvers::report
        (
            "smoothSolver " + fieldName_
          + ": fixed-iteration path, residual after the last sweep only"
        );
    }

    // As smoothSolver::solve with nSweeps_ >= 0, minus the residual norm
    // after every sweep but the last (checkConvergence is always false).
    PrecisionAdaptor<solveScalar, scalar> tpsi(psi_s);
    solveScalarField& psi = tpsi.ref();

    solverPerformance solverPerf(typeName, fieldName_);

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
