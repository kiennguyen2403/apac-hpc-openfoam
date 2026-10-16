/*---------------------------------------------------------------------------*\
    Replacement for linearViscousStress<incompressible RAS>::divDevRhoReff(U)
    (explicit specialisation; reached through the vtable under LD_PRELOAD).

    Stock:
        - fvc::div((alpha*rho*nuEff())*dev2(T(fvc::grad(U))))
        - fvm::laplacian(alpha*rho*nuEff(), U)

    The explicit term makes three cell tensor fields (T, dev2, nuEff*...),
    a face vector field (Sf & linear interpolate) and then integrates it over
    the cells (surfaceIntegrate), and nuEff() is evaluated twice. Here:
        - nuEff() once,
        - X = nuEff*dev2(T(gradU)) in one cell loop (same per-cell functions),
        - one face loop: Sf & (lambda*(X[P] - X[N]) + X[N]) added to the owner
          and subtracted from the neighbour, in face order, then the patches
          in patch order, then divided by the cell volumes: the order of the
          stock surfaceIntegrate, so every cell sum is unchanged.
    The implicit term and the combination are the stock code. Bit-identical.

    Only for div scheme "Gauss linear" and patches that are either not coupled
    or processor patches (whose neighbour values are the stored patch values);
    otherwise, or with STEP_REUSE_OFF=divdev, the stock function runs.
\*---------------------------------------------------------------------------*/

#include "IncompressibleTurbulenceModel.H"
#include "incompressible/transportModel/transportModel.H"
#include "RASModel.H"
#include "linearViscousStress.H"
#include "fvc.H"
#include "fvm.H"
#include "processorFvPatch.H"
#include "stepReuse.H"

namespace Foam
{

typedef RASModel<IncompressibleTurbulenceModel<transportModel>>
    stepReuseIncompressibleRAS;

template<>
tmp<fvVectorMatrix>
linearViscousStress<stepReuseIncompressibleRAS>::divDevRhoReff
(
    volVectorField& U
) const
{
    const fvMesh& mesh = this->mesh_;

    tmp<volScalarField> tnuEff(this->nuEff());
    const volScalarField& nuEff = tnuEff();

    // Name of the stock explicit field, (nuEff*dev2(T(grad(U)))), and its
    // div scheme
    const word Xname
    (
        "(" + nuEff.name() + "*dev2(T(grad(" + U.name() + "))))"
    );
    const word divName("div(" + Xname + ')');

    bool fused = !stepReuse::off("divdev");

    if (fused)
    {
        const ITstream& is = mesh.divScheme(divName);
        fused =
            is.size() == 2
         && is[0].isWord() && is[0].wordToken() == "Gauss"
         && is[1].isWord() && is[1].wordToken() == "linear";
    }

    if (fused)
    {
        for (const fvPatch& p : mesh.boundary())
        {
            if (p.coupled() && !isA<processorFvPatch>(p))
            {
                fused = false;
            }
        }
    }

    if (!fused)
    {
        stepReuse::report("divdev", "divDevReff: stock");

        typedef tmp<fvVectorMatrix> (*divDevFn)
        (
            const linearViscousStress<stepReuseIncompressibleRAS>*,
            volVectorField&
        );
        static const divDevFn stock = reinterpret_cast<divDevFn>
        (
            stepReuse::next
            (
                "_ZNK4Foam19linearViscousStressINS_8RASModelINS_"
                "29IncompressibleTurbulenceModelINS_14transportModelEEEEEE"
                "13divDevRhoReffERNS_14GeometricFieldINS_6VectorIdEENS_"
                "12fvPatchFieldENS_7volMeshEEE"
            )
        );
        return stock(this, U);
    }

    stepReuse::report("divdev", "divDevReff: fused explicit term");

    tmp<volTensorField> tgradU(fvc::grad(U));
    const volTensorField& gradU = tgradU();

    const scalarField& nu = nuEff.primitiveField();
    const tensorField& G = gradU.primitiveField();

    // X = nuEff*dev2(T(gradU)), cell values
    tensorField X(G.size());
    forAll(X, celli)
    {
        X[celli] = nu[celli]*dev2(T(G[celli]));
    }

    // Explicit term, as -fvc::div(X) with Gauss linear
    const surfaceScalarField& lambdas = mesh.surfaceInterpolation::weights();
    const surfaceVectorField& Sf = mesh.Sf();

    const scalarField& lambda = lambdas.primitiveField();
    const vectorField& Sfi = Sf.primitiveField();
    const labelUList& P = mesh.owner();
    const labelUList& N = mesh.neighbour();

    tmp<volVectorField> tnegDiv
    (
        new volVectorField
        (
            IOobject
            (
                divName,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh,
            dimensioned<vector>
            (
                Sf.dimensions()*nuEff.dimensions()*gradU.dimensions()/dimVol,
                Zero
            )
        )
    );
    vectorField& ivf = tnegDiv.ref().primitiveFieldRef();

    forAll(P, facei)
    {
        const vector sfi =
            Sfi[facei] & (lambda[facei]*(X[P[facei]] - X[N[facei]]) + X[N[facei]]);

        ivf[P[facei]] += sfi;
        ivf[N[facei]] -= sfi;
    }

    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& p = mesh.boundary()[patchi];
        const labelUList& faceCells = p.faceCells();

        // Patch values of X as the stock field operations make them
        const tensorField Xb
        (
            nuEff.boundaryField()[patchi]
           *dev2(gradU.boundaryField()[patchi].T())
        );

        vectorField psf;

        if (p.coupled())
        {
            // processor: neighbour values are the patch values
            tensorField Xi(faceCells.size());
            forAll(faceCells, facei)
            {
                Xi[facei] = X[faceCells[facei]];
            }

            psf =
                Sf.boundaryField()[patchi]
              & lerp(Xb, Xi, lambdas.boundaryField()[patchi]);
        }
        else
        {
            psf = Sf.boundaryField()[patchi] & Xb;
        }

        forAll(faceCells, facei)
        {
            ivf[faceCells[facei]] += psf[facei];
        }
    }

    ivf /= mesh.Vsc();
    ivf.negate();

    X.clear();
    tgradU.clear();

    return tnegDiv - fvm::laplacian(nuEff, U);
}

} // End namespace Foam


// ************************************************************************* //
