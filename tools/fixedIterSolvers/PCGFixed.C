#include "PCGFixed.H"
#include "fixedIterControl.H"
#include "PrecisionAdaptor.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(PCGFixed, 0);

    lduMatrix::solver::addsymMatrixConstructorToTable<PCGFixed>
        addPCGFixedSymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PCGFixed::PCGFixed
(
    const word& fieldName,
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const FieldField<Field, scalar>& interfaceIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const dictionary& solverControls
)
:
    lduMatrix::solver
    (
        fieldName,
        matrix,
        interfaceBouCoeffs,
        interfaceIntCoeffs,
        interfaces,
        solverControls
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::solverPerformance Foam::PCGFixed::scalarSolve
(
    solveScalarField& psi,
    const solveScalarField& source,
    const direction cmpt
) const
{
    // --- Setup class containing solver performance data
    solverPerformance solverPerf
    (
        lduMatrix::preconditioner::getName(controlDict_) + word("PCG"),
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

    if ((log_ >= 2) || (lduMatrix::debug >= 2))
    {
        Info<< "   Normalisation factor = " << normFactor << endl;
    }

    // --- Calculate normalised residual norm
    solverPerf.initialResidual() =
        gSumMag(rA, matrix().mesh().comm())
       /normFactor;
    solverPerf.finalResidual() = solverPerf.initialResidual();

    // --- Check convergence, solve if not converged
    if
    (
        minIter_ > 0
     || !solverPerf.checkConvergence(tolerance_, relTol_, log_)
    )
    {
        const bool fixedIter =
            convergenceImpossible(tolerance_, relTol_, log_);

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

            if (!fixedIter)
            {
                solverPerf.finalResidual() =
                    gSumMag(rA, matrix().mesh().comm())
                   /normFactor;
            }

        } while
        (
            fixedIter
          ? (
                ++solverPerf.nIterations() < maxIter_
             || solverPerf.nIterations() < minIter_
            )
          : (
                (
                  ++solverPerf.nIterations() < maxIter_
                && !solverPerf.checkConvergence(tolerance_, relTol_, log_)
                )
             || solverPerf.nIterations() < minIter_
            )
        );

        // rA is unchanged since the last update (also after a singularity
        // break), so this equals PCG's last per-iteration norm
        if (fixedIter)
        {
            solverPerf.finalResidual() =
                gSumMag(rA, matrix().mesh().comm())
               /normFactor;
        }
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


Foam::solverPerformance Foam::PCGFixed::solve
(
    scalarField& psi_s,
    const scalarField& source,
    const direction cmpt
) const
{
    PrecisionAdaptor<solveScalar, scalar> tpsi(psi_s);
    return scalarSolve
    (
        tpsi.ref(),
        ConstPrecisionAdaptor<solveScalar, scalar>(source)(),
        cmpt
    );
}


// ************************************************************************* //
