/*---------------------------------------------------------------------------*\
    LD_PRELOAD replacement for Foam::PCG::scalarSolve, the coarsest-level
    solver of the shipped GAMG settings for p:

        coarsestLevelCorr { solver PCG; preconditioner DIC;
                            minIter 35; maxIter 35; }

    With minIter == maxIter the iteration count is fixed whatever
    checkConvergence returns (the loop runs while nIterations < minIter), and
    with tolerance 0 / relTol 0 convergence is impossible. In both cases the
    residual norm after every iteration (a global sum over all ranks) only
    feeds a test that cannot change the count. Here it is computed once,
    after the last iteration: same iterations, same psi, same initial and
    final residuals, same singularity handling. At 4 nodes: 12 V-cycles x 34
    norms = 408 fewer Allreduces per time step (the census, job 182140240,
    counted 3 global sums per coarsest-level iteration).

    Any other controls, per-iteration printing (log >= 2 or debug >= 2), or
    STEP_REUSE_OFF=pcg run the stock function (next definition, in
    libOpenFOAM). GAMG only prints the coarsest solve's performance in debug
    mode, so the converged flag this path leaves unset is never read.
\*---------------------------------------------------------------------------*/

#include "PCG.H"
#include "PrecisionAdaptor.H"
#include "fixedIterControl.H"
#include "fixedIterSolvers.H"


Foam::solverPerformance Foam::PCG::scalarSolve
(
    solveScalarField& psi,
    const solveScalarField& source,
    const direction cmpt
) const
{
    const bool fixedCount =
    (
        minIter_ > 0
     && maxIter_ <= minIter_
     && log_ < 2
     && lduMatrix::debug < 2
     && solverPerformance::debug < 2
    )
     || convergenceImpossible(tolerance_, relTol_, log_);

    if (!fixedCount || fixedIterSolvers::off("pcg"))
    {
        // Itanium C++ ABI: the hidden return-value address comes first, then
        // this, as for a free function taking the object as first argument
        typedef solverPerformance (*solveFn)
        (
            const PCG*,
            solveScalarField&,
            const solveScalarField&,
            const direction
        );
        static const solveFn stock = reinterpret_cast<solveFn>
        (
            fixedIterSolvers::next
            (
                "_ZNK4Foam3PCG11scalarSolveERNS_5FieldIdEERKS2_h"
            )
        );

        if (fixedIterSolvers::reporting())
        {
            fixedIterSolvers::report("PCG " + fieldName_ + ": stock solve");
        }

        return stock(this, psi, source, cmpt);
    }

    if (fixedIterSolvers::reporting())
    {
        fixedIterSolvers::report
        (
            "PCG " + fieldName_ + ": fixed iteration count (minIter "
          + Foam::name(minIter_) + ", maxIter " + Foam::name(maxIter_)
          + "), residual norm after the last iteration only"
        );
    }

    // As PCG::scalarSolve, minus the residual norm inside the loop
    solverPerformance solverPerf
    (
        lduMatrix::preconditioner::getName(controlDict_) + typeName,
        fieldName_
    );

    label nCells = psi.size();

    solveScalar* __restrict__ psiPtr = psi.begin();

    solveScalarField pA(nCells);
    solveScalar* __restrict__ pAPtr = pA.begin();

    solveScalarField wA(nCells);
    solveScalar* __restrict__ wAPtr = wA.begin();

    solveScalar wArA = solverPerf.great_;
    solveScalar wArAold = wArA;

    // --- Calculate A.psi
    matrix_.Amul(wA, psi, interfaceBouCoeffs_, interfaces_, cmpt);

    // --- Calculate initial residual field
    solveScalarField rA(source - wA);
    solveScalar* __restrict__ rAPtr = rA.begin();

    matrix().setResidualField
    (
        ConstPrecisionAdaptor<scalar, solveScalar>(rA)(),
        fieldName_,
        true
    );

    // --- Calculate normalisation factor
    solveScalar normFactor = this->normFactor(psi, source, wA, pA);

    // --- Calculate normalised residual norm
    solverPerf.initialResidual() =
        gSumMag(rA, matrix().mesh().comm())
       /normFactor;
    solverPerf.finalResidual() = solverPerf.initialResidual();

    // --- The count is fixed: solve (stock: minIter_ > 0 || !converged)
    {
        // --- Select and construct the preconditioner
        if (!preconPtr_)
        {
            preconPtr_ = lduMatrix::preconditioner::New
            (
                *this,
                controlDict_
            );
        }

        // --- Solver iteration
        do
        {
            // --- Store previous wArA
            wArAold = wArA;

            // --- Precondition residual
            preconPtr_->precondition(wA, rA, cmpt);

            // --- Update search directions:
            wArA = gSumProd(wA, rA, matrix().mesh().comm());

            if (solverPerf.nIterations() == 0)
            {
                for (label cell=0; cell<nCells; cell++)
                {
                    pAPtr[cell] = wAPtr[cell];
                }
            }
            else
            {
                const solveScalar beta = wArA/wArAold;

                for (label cell=0; cell<nCells; cell++)
                {
                    pAPtr[cell] = wAPtr[cell] + beta*pAPtr[cell];
                }
            }


            // --- Update preconditioned residual
            matrix_.Amul(wA, pA, interfaceBouCoeffs_, interfaces_, cmpt);

            solveScalar wApA = gSumProd(wA, pA, matrix().mesh().comm());

            // --- Test for singularity
            if (solverPerf.checkSingularity(mag(wApA)/normFactor)) break;


            // --- Update solution and residual:

            const solveScalar alpha = wArA/wApA;

            for (label cell=0; cell<nCells; cell++)
            {
                psiPtr[cell] += alpha*pAPtr[cell];
                rAPtr[cell] -= alpha*wAPtr[cell];
            }

            // (stock: finalResidual = gSumMag(rA)/normFactor here, then
            //  checkConvergence in the loop condition; both without effect
            //  on the count)
        } while
        (
            ++solverPerf.nIterations() < maxIter_
         || solverPerf.nIterations() < minIter_
        );

        // --- The norm the stock loop computed after its last update (the
        //     same reduction over the same rA: bit-identical)
        solverPerf.finalResidual() =
            gSumMag(rA, matrix().mesh().comm())
           /normFactor;
    }

    if (preconPtr_)
    {
        preconPtr_->setFinished(solverPerf);
    }

    matrix().setResidualField
    (
        ConstPrecisionAdaptor<scalar, solveScalar>(rA)(),
        fieldName_,
        false
    );

    return solverPerf;
}


// ************************************************************************* //
