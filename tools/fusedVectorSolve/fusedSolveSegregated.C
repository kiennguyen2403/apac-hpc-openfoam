/*---------------------------------------------------------------------------*\
    fusedSolveSegregated: replacement for fvMatrix<vector>::solveSegregated,
    loaded with LD_PRELOAD (libfiniteVolume exports the generic instantiation
    as a weak symbol called through the PLT, so this definition takes over).

    For a segregated vector solve with smoothSolver/smoothSolverFixed +
    GaussSeidel and fixed-iteration controls (tolerance 0, relTol 0), the
    three components are independent and share the matrix off-diagonals and
    addressing; only the boundary-corrected diagonal, the source and the
    interface coefficients differ. Each component is set up exactly as in the
    generic code (own diagonal copy, source, interface correction, initial
    residual and norm factor), then one GaussSeidel sweep updates all three
    components, reading upper/lower/addressing once instead of three times.
    Per-component arithmetic and order are unchanged, so results, residuals
    and log output are bit-identical. Any other case runs the generic code
    (copied verbatim).
\*---------------------------------------------------------------------------*/

#include "fvMatrices.H"
#include "volFields.H"
#include "PrecisionAdaptor.H"
#include "FixedList.H"
#include "PtrList.H"

#include <cstdlib>
#include <iostream>

namespace Foam
{
namespace fusedVectorSolve
{

struct controls
{
    label nSweeps = 1;
    label maxIter = lduMatrix::defaultMaxIter;
    label minIter = 0;
    bool ok = false;
};

// Controls for which smoothSolver + GaussSeidel can never converge early and
// print nothing per iteration, so the fused sweep reproduces it exactly
static controls readControls(const dictionary& d)
{
    controls c;
    word solver;
    if (!d.readIfPresent("solver", solver)) return c;
    if (solver != "smoothSolver" && solver != "smoothSolverFixed") return c;
    if (!d.found("smoother", keyType::LITERAL)) return c;
    if (lduMatrix::smoother::getName(d) != "GaussSeidel") return c;
    if (d.found("norm", keyType::LITERAL)) return c;

    int log = 1;
    scalar tolerance = lduMatrix::defaultTolerance;
    scalar relTol = 0;
    d.readIfPresent("log", log);
    d.readIfPresent("minIter", c.minIter);
    d.readIfPresent("maxIter", c.maxIter);
    d.readIfPresent("tolerance", tolerance);
    d.readIfPresent("relTol", relTol);
    c.nSweeps = d.getOrDefault<label>("nSweeps", 1);

    c.ok =
        c.nSweeps > 0
     && tolerance <= 0
     && relTol <= solverPerformance::small_
     && log < 2
     && lduMatrix::debug < 2
     && solverPerformance::debug < 2;

    return c;
}

} // End namespace fusedVectorSolve
} // End namespace Foam


template<>
Foam::SolverPerformance<Foam::vector>
Foam::fvMatrix<Foam::vector>::solveSegregated
(
    const dictionary& solverControls
)
{
    typedef vector Type;

    if (useImplicit_)
    {
        FatalErrorInFunction
            << "Implicit option is not allowed for type: " << Type::typeName
            << exit(FatalError);
    }

    if (debug)
    {
        Info.masterStream(this->mesh().comm())
            << "fvMatrix<Type>::solveSegregated"
               "(const dictionary& solverControls) : "
               "solving fvMatrix<Type>"
            << endl;
    }

    const int logLevel =
        solverControls.getOrDefault<int>
        (
            "log",
            SolverPerformance<Type>::debug
        );

    auto& psi = psi_.constCast();

    SolverPerformance<Type> solverPerfVec
    (
        "fvMatrix<Type>::solveSegregated",
        psi.name()
    );

    scalarField saveDiag(diag());

    Field<Type> source(source_);

    // At this point include the boundary source from the coupled boundaries.
    // This is corrected for the implicit part by updateMatrixInterfaces within
    // the component loop.
    addBoundarySource(source);

    typename Type::labelType validComponents
    (
        psi.mesh().template validComponents<Type>()
    );

    const fusedVectorSolve::controls ctl =
        fusedVectorSolve::readControls(solverControls);

    bool fused = ctl.ok;
    for (direction cmpt=0; cmpt<Type::nComponents; cmpt++)
    {
        if (validComponents[cmpt] == -1) fused = false;
    }

    // Opt-in check that this library is active (stderr, once per process):
    // FUSED_VECTOR_SOLVE_REPORT=1
    static bool reported = false;
    if (!reported && std::getenv("FUSED_VECTOR_SOLVE_REPORT") && UPstream::master())
    {
        reported = true;
        std::cerr << "fusedVectorSolve: " << psi.name()
            << (fused ? " fused path" : " generic path") << std::endl;
    }

    if (fused)
    {
        const label comm = this->mesh().comm();
        const label nCells = psi.primitiveField().size();

        lduInterfaceFieldPtrsList interfaces =
            psi.boundaryField().scalarInterfaces();

        FixedList<scalarField, 3> psiC;
        FixedList<scalarField, 3> diagC;
        FixedList<scalarField, 3> sourceC;
        PtrList<FieldField<Field, scalar>> bouC(3);
        FixedList<solverPerformance, 3> perf;
        FixedList<scalar, 3> normFactor;

        // --- per-component setup and initial residual, as the generic loop
        //     followed by smoothSolver::solve up to the smoothing loop
        for (direction cmpt=0; cmpt<3; cmpt++)
        {
            const word fieldName(psi.name() + pTraits<Type>::componentNames[cmpt]);

            psiC[cmpt] = psi.primitiveField().component(cmpt);
            diagC[cmpt] = saveDiag;
            addBoundaryDiag(diagC[cmpt], cmpt);
            sourceC[cmpt] = source.component(cmpt);
            bouC.set(cmpt, new FieldField<Field, scalar>(boundaryCoeffs_.component(cmpt)));

            {
                PrecisionAdaptor<solveScalar, scalar> sourceCmpt_ss(sourceC[cmpt]);
                ConstPrecisionAdaptor<solveScalar, scalar> psiCmpt_ss(psiC[cmpt]);

                const label startRequest = UPstream::nRequests();

                initMatrixInterfaces
                (
                    true,
                    bouC[cmpt],
                    interfaces,
                    psiCmpt_ss(),
                    sourceCmpt_ss.ref(),
                    cmpt
                );

                updateMatrixInterfaces
                (
                    true,
                    bouC[cmpt],
                    interfaces,
                    psiCmpt_ss(),
                    sourceCmpt_ss.ref(),
                    cmpt,
                    startRequest
                );
            }

            diag() = diagC[cmpt];

            perf[cmpt] = solverPerformance(word("smoothSolver"), fieldName);

            solveScalarField Apsi(nCells);
            solveScalarField temp(nCells);

            this->lduMatrix::Amul(Apsi, psiC[cmpt], bouC[cmpt], interfaces, cmpt);

            // lduMatrix::solver::normFactor (DEFAULT_NORM)
            this->lduMatrix::sumA(temp, bouC[cmpt], interfaces);
            temp *= gAverage(psiC[cmpt], comm);
            normFactor[cmpt] =
                gSum
                (
                    (mag(Apsi - temp) + mag(sourceC[cmpt] - temp))(),
                    comm
                ) + solverPerformance::small_;

            solveScalarField residual(sourceC[cmpt] - Apsi);

            this->lduMatrix::setResidualField(residual, fieldName, true);

            perf[cmpt].initialResidual() =
                gSumMag(residual, comm)/normFactor[cmpt];
            perf[cmpt].finalResidual() = perf[cmpt].initialResidual();
        }

        // --- fused GaussSeidel sweeps (GaussSeidelSmoother::smooth for the
        //     three components in one pass over the matrix)
        {
            // const access: the non-const lower() would allocate a separate
            // lower array for a symmetric matrix
            const lduMatrix& constMatrix = *this;
            const scalar* const __restrict__ upperPtr = constMatrix.upper().begin();
            const scalar* const __restrict__ lowerPtr = constMatrix.lower().begin();
            const label* const __restrict__ uPtr = lduAddr().upperAddr().begin();
            const label* const __restrict__ ownStartPtr =
                lduAddr().ownerStartAddr().begin();

            FixedList<solveScalarField, 3> bPrime;
            for (direction cmpt=0; cmpt<3; cmpt++) bPrime[cmpt].resize(nCells);

            solveScalar* __restrict__ psi0 = psiC[0].begin();
            solveScalar* __restrict__ psi1 = psiC[1].begin();
            solveScalar* __restrict__ psi2 = psiC[2].begin();
            solveScalar* __restrict__ b0 = bPrime[0].begin();
            solveScalar* __restrict__ b1 = bPrime[1].begin();
            solveScalar* __restrict__ b2 = bPrime[2].begin();
            const scalar* const __restrict__ d0 = diagC[0].begin();
            const scalar* const __restrict__ d1 = diagC[1].begin();
            const scalar* const __restrict__ d2 = diagC[2].begin();

            label nIterations = 0;
            do
            {
                for (label sweep=0; sweep<ctl.nSweeps; sweep++)
                {
                    for (direction cmpt=0; cmpt<3; cmpt++)
                    {
                        bPrime[cmpt] = sourceC[cmpt];

                        const label startRequest = UPstream::nRequests();

                        initMatrixInterfaces
                        (
                            false,
                            bouC[cmpt],
                            interfaces,
                            psiC[cmpt],
                            bPrime[cmpt],
                            cmpt
                        );

                        updateMatrixInterfaces
                        (
                            false,
                            bouC[cmpt],
                            interfaces,
                            psiC[cmpt],
                            bPrime[cmpt],
                            cmpt,
                            startRequest
                        );
                    }

                    solveScalar p0, p1, p2;
                    label fStart;
                    label fEnd = ownStartPtr[0];

                    for (label celli=0; celli<nCells; celli++)
                    {
                        fStart = fEnd;
                        fEnd = ownStartPtr[celli + 1];

                        p0 = b0[celli];
                        p1 = b1[celli];
                        p2 = b2[celli];

                        for (label facei=fStart; facei<fEnd; facei++)
                        {
                            const label j = uPtr[facei];
                            p0 -= upperPtr[facei]*psi0[j];
                            p1 -= upperPtr[facei]*psi1[j];
                            p2 -= upperPtr[facei]*psi2[j];
                        }

                        p0 /= d0[celli];
                        p1 /= d1[celli];
                        p2 /= d2[celli];

                        for (label facei=fStart; facei<fEnd; facei++)
                        {
                            const label j = uPtr[facei];
                            b0[j] -= lowerPtr[facei]*p0;
                            b1[j] -= lowerPtr[facei]*p1;
                            b2[j] -= lowerPtr[facei]*p2;
                        }

                        psi0[celli] = p0;
                        psi1[celli] = p1;
                        psi2[celli] = p2;
                    }
                }
            } while
            (
                (nIterations += ctl.nSweeps) < ctl.maxIter
             || nIterations < ctl.minIter
            );

            for (direction cmpt=0; cmpt<3; cmpt++)
            {
                perf[cmpt].nIterations() = nIterations;
            }
        }

        // --- final residuals, as smoothSolver::solve after the loop
        for (direction cmpt=0; cmpt<3; cmpt++)
        {
            diag() = diagC[cmpt];

            solveScalarField residual
            (
                this->lduMatrix::residual
                (
                    psiC[cmpt],
                    sourceC[cmpt],
                    bouC[cmpt],
                    interfaces,
                    cmpt
                )
            );

            perf[cmpt].finalResidual() =
                gSumMag(residual, comm)/normFactor[cmpt];

            this->lduMatrix::setResidualField
            (
                residual,
                psi.name() + pTraits<Type>::componentNames[cmpt],
                false
            );
        }

        // --- report and store in component order, as the generic loop
        for (direction cmpt=0; cmpt<3; cmpt++)
        {
            if (logLevel)
            {
                perf[cmpt].print(Info.masterStream(this->mesh().comm()));
            }

            solverPerfVec.replace(cmpt, perf[cmpt]);
            solverPerfVec.solverName() = perf[cmpt].solverName();

            psi.primitiveFieldRef().replace(cmpt, psiC[cmpt]);
        }
        diag() = saveDiag;
    }
    else
    {
        // --- generic code (fvMatrixSolve.C, v2512), unchanged
        for (direction cmpt=0; cmpt<Type::nComponents; cmpt++)
        {
            if (validComponents[cmpt] == -1) continue;

            // copy field and source

            scalarField psiCmpt(psi.primitiveField().component(cmpt));
            addBoundaryDiag(diag(), cmpt);

            scalarField sourceCmpt(source.component(cmpt));

            FieldField<Field, scalar> bouCoeffsCmpt
            (
                boundaryCoeffs_.component(cmpt)
            );

            FieldField<Field, scalar> intCoeffsCmpt
            (
                internalCoeffs_.component(cmpt)
            );

            lduInterfaceFieldPtrsList interfaces =
                psi.boundaryField().scalarInterfaces();

            // Use the initMatrixInterfaces and updateMatrixInterfaces to correct
            // bouCoeffsCmpt for the explicit part of the coupled boundary
            // conditions
            {
                PrecisionAdaptor<solveScalar, scalar> sourceCmpt_ss(sourceCmpt);
                ConstPrecisionAdaptor<solveScalar, scalar> psiCmpt_ss(psiCmpt);

                const label startRequest = UPstream::nRequests();

                initMatrixInterfaces
                (
                    true,
                    bouCoeffsCmpt,
                    interfaces,
                    psiCmpt_ss(),
                    sourceCmpt_ss.ref(),
                    cmpt
                );

                updateMatrixInterfaces
                (
                    true,
                    bouCoeffsCmpt,
                    interfaces,
                    psiCmpt_ss(),
                    sourceCmpt_ss.ref(),
                    cmpt,
                    startRequest
                );
            }

            solverPerformance solverPerf;

            // Solver call
            solverPerf = lduMatrix::solver::New
            (
                psi.name() + pTraits<Type>::componentNames[cmpt],
                *this,
                bouCoeffsCmpt,
                intCoeffsCmpt,
                interfaces,
                solverControls
            )->solve(psiCmpt, sourceCmpt, cmpt);

            if (logLevel)
            {
                solverPerf.print(Info.masterStream(this->mesh().comm()));
            }

            solverPerfVec.replace(cmpt, solverPerf);
            solverPerfVec.solverName() = solverPerf.solverName();

            psi.primitiveFieldRef().replace(cmpt, psiCmpt);
            diag() = saveDiag;
        }
    }

    psi.correctBoundaryConditions();

    psi.mesh().data().setSolverPerformance(psi.name(), solverPerfVec);

    return solverPerfVec;
}
