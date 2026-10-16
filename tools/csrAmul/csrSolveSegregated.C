/*---------------------------------------------------------------------------*\
    csrAmul: replacement for fvMatrix<scalar>::solveSegregated (fvScalarMatrix.C,
    v2512), loaded with LD_PRELOAD. Identical to the original except that it
    builds lduMatrix::lowerCSR() before the solver is called, which makes
    lduMatrix::Amul use its built-in cell-based loop (gather per cell instead
    of scatter per face) for the finest-level products (e.g. GAMG residuals and
    correction scaling of p). Same arithmetic in the same order (losort is
    stable), so results are bit-identical.
\*---------------------------------------------------------------------------*/

#include "fvScalarMatrix.H"
#include "extrapolatedCalculatedFvPatchFields.H"
#include "profiling.H"
#include "PrecisionAdaptor.H"
#include "jumpCyclicFvPatchField.H"
#include "cyclicPolyPatch.H"
#include "cyclicAMIPolyPatch.H"

#include <cstdlib>
#include <iostream>

template<>
Foam::solverPerformance Foam::fvMatrix<Foam::scalar>::solveSegregated
(
    const dictionary& solverControls
)
{
    if (debug)
    {
        Info.masterStream(this->mesh().comm())
            << "fvMatrix<scalar>::solveSegregated"
               "(const dictionary& solverControls) : "
               "solving fvMatrix<scalar>"
            << endl;
    }

    const int logLevel =
        solverControls.getOrDefault<int>
        (
            "log",
            solverPerformance::debug
        );

    scalarField saveLower;
    scalarField saveUpper;

    if (useImplicit_)
    {
        createOrUpdateLduPrimitiveAssembly();

        if (psi_.mesh().fluxRequired(psi_.name()))
        {
            // Save lower/upper for flux calculation
            if (asymmetric())
            {
                saveLower = lower();
            }
            saveUpper = upper();
        }

        setLduMesh(*lduMeshPtr());
        transferFvMatrixCoeffs();
        setBounAndInterCoeffs();
        direction cmpt = 0;
        manipulateMatrix(cmpt);
    }

    scalarField saveDiag(diag());
    addBoundaryDiag(diag(), 0);

    scalarField totalSource(source_);
    addBoundarySource(totalSource, false);

    lduInterfaceFieldPtrsList interfaces;
    PtrDynList<lduInterfaceField> newInterfaces;
    if (!useImplicit_)
    {
        interfaces = this->psi(0).boundaryField().scalarInterfaces();
    }
    else
    {
        setInterfaces(interfaces, newInterfaces);
    }

    tmp<scalarField> tpsi;
    if (!useImplicit_)
    {
        tpsi.ref(psi_.constCast().primitiveFieldRef());
    }
    else
    {
        tpsi = tmp<scalarField>::New(lduAddr().size(), Zero);
        scalarField& psi = tpsi.ref();

        for (label fieldi = 0; fieldi < nMatrices(); fieldi++)
        {
            const label cellOffset = lduMeshPtr()->cellOffsets()[fieldi];
            const auto& psiInternal = this->psi(fieldi).primitiveField();

            forAll(psiInternal, localCellI)
            {
                psi[cellOffset + localCellI] = psiInternal[localCellI];
            }
        }
    }
    scalarField& psi = tpsi.ref();

    // csrAmul: build the cell-ordered copy of the lower coefficients so that
    // lduMatrix::Amul takes its cell-based (gather) path for this matrix.
    // losort is stable, so the gather visits each cell's contributions in the
    // same order as the face loop: bit-identical results.
    if (!useImplicit_)
    {
        (void)static_cast<const lduMatrix&>(*this).lowerCSR();
    }
    const bool csrBefore = hasLowerCSR();

    // Solver call
    solverPerformance solverPerf = lduMatrix::solver::New
    (
        this->psi(0).name(),
        *this,
        boundaryCoeffs_,
        internalCoeffs_,
        interfaces,
        solverControls
    )->solve(psi, totalSource);

    if (useImplicit_)
    {
        for (label fieldi = 0; fieldi < nMatrices(); fieldi++)
        {
            auto& psiInternal =
                this->psi(fieldi).constCast().primitiveFieldRef();

            const label cellOffset = lduMeshPtr()->cellOffsets()[fieldi];

            forAll(psiInternal, localCellI)
            {
                psiInternal[localCellI] = psi[localCellI + cellOffset];
            }
        }
    }

    // Opt-in check (stderr, once per process): CSR_AMUL_REPORT=1. "kept" means
    // nothing reset the CSR copy during the solve (non-const upper()/lower()).
    static bool reported = false;
    if (!reported && std::getenv("CSR_AMUL_REPORT") && UPstream::master())
    {
        reported = true;
        std::cerr << "csrAmul: " << this->psi(0).name()
            << " lowerCSR " << (csrBefore ? "built" : "not built")
            << ", after solve " << (hasLowerCSR() ? "kept" : "reset") << std::endl;
    }

    if (logLevel)
    {
        solverPerf.print(Info.masterStream(mesh().comm()));
    }

    diag() = saveDiag;

    if (useImplicit_)
    {
        if (psi_.mesh().fluxRequired(psi_.name()))
        {
            // Restore lower/upper
            if (asymmetric())
            {
                lower().setSize(saveLower.size());
                lower() = saveLower;
            }

            upper().setSize(saveUpper.size());
            upper() = saveUpper;
        }
        // Set the original lduMesh
        setLduMesh(psi_.mesh());
    }

    for (label fieldi = 0; fieldi < nMatrices(); fieldi++)
    {
        auto& localPsi = this->psi(fieldi).constCast();

        localPsi.correctBoundaryConditions();
        localPsi.mesh().data().setSolverPerformance
        (
            localPsi.name(),
            solverPerf
        );
    }

    return solverPerf;
}
