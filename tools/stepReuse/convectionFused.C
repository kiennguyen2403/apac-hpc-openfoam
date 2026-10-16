/*---------------------------------------------------------------------------*\
    Replacement for fv::gaussConvectionScheme<vector>::fvmDiv with the
    linearUpwindV interpolation (div(phi,U) Gauss linearUpwindV grad(U),
    also under 'bounded', which wraps it).

    Stock: weights = pos0(phi) (a face field), lower = -weights*phi, upper =
    lower + phi, negSumDiag, patch coefficients; then the explicit
    correction: linearUpwindV::correction(U) (a face vector field), phi times
    it (a face vector field), surfaceIntegrate (a cell vector field), source
    -= V*that.

    Here one loop over the internal faces computes the weight, lower, upper,
    the diagonal (owner then neighbour, in face order, as negSumDiag) and the
    limited correction (the linearUpwindV code, unchanged), and adds phi
    times it to the owner and subtracts it from the neighbour (face order,
    then the patches, as surfaceIntegrate). Same expressions on the same
    values: bit-identical.

    The gradient is the one linearUpwindV uses (its grad scheme name, read
    from the div scheme entry, through the solution cache). Only for
    "[bounded] Gauss linearUpwindV <grad>"; otherwise, or with
    STEP_REUSE_OFF=conv, the stock function runs.
\*---------------------------------------------------------------------------*/

#include "gaussConvectionScheme.H"
#include "gradScheme.H"
#include "fvMatrices.H"
#include "stepReuse.H"

namespace Foam
{
namespace stepReuseDetail
{

// Grad scheme name in "[bounded] Gauss linearUpwindV <grad>", else empty
static word linearUpwindVGradName(const ITstream& is)
{
    label i = 0;
    if (is.size() == 4 && is[0].isWord() && is[0].wordToken() == "bounded")
    {
        i = 1;
    }
    if
    (
        is.size() != 3 + i
     || !(is[i].isWord() && is[i].wordToken() == "Gauss")
     || !(is[i+1].isWord() && is[i+1].wordToken() == "linearUpwindV")
     || !is[i+2].isWord()
    )
    {
        return word::null;
    }
    return is[i+2].wordToken();
}

} // End namespace stepReuseDetail


template<>
tmp<fvMatrix<vector>> fv::gaussConvectionScheme<vector>::fvmDiv
(
    const surfaceScalarField& faceFlux,
    const volVectorField& vf
) const
{
    const fvMesh& mesh = this->mesh();

    word gradName;
    if (!stepReuse::off("conv") && tinterpScheme_().type() == "linearUpwindV")
    {
        gradName = stepReuseDetail::linearUpwindVGradName
        (
            mesh.divScheme("div(" + faceFlux.name() + ',' + vf.name() + ')')
        );
    }

    if (gradName.empty())
    {
        typedef tmp<fvMatrix<vector>> (*fn)
        (
            const fv::gaussConvectionScheme<vector>*,
            const surfaceScalarField&,
            const volVectorField&
        );
        static const fn stock = reinterpret_cast<fn>
        (
            stepReuse::next
            (
                "_ZNK4Foam2fv21gaussConvectionSchemeINS_6VectorIdEEE6fvmDivERKNS_"
                "14GeometricFieldIdNS_13fvsPatchFieldENS_11surfaceMeshEEERKNS5_"
                "IS3_NS_12fvPatchFieldENS_7volMeshEEE"
            )
        );
        return stock(this, faceFlux, vf);
    }

    stepReuse::report("conv", "convection: fused Gauss linearUpwindV");

    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();

    tmp<fvMatrix<vector>> tfvm
    (
        new fvMatrix<vector>
        (
            vf,
            faceFlux.dimensions()*vf.dimensions()
        )
    );
    fvMatrix<vector>& fvm = tfvm.ref();

    scalarField& lower = fvm.lower();
    scalarField& upper = fvm.upper();
    scalarField& diag = fvm.diag();

    // Patch coefficients, with the upwind weights pos0(phi) of the patch
    forAll(vf.boundaryField(), patchi)
    {
        const fvPatchField<vector>& psf = vf.boundaryField()[patchi];
        const fvsPatchScalarField& patchFlux = faceFlux.boundaryField()[patchi];
        const scalarField pw(pos0(patchFlux));

        fvm.internalCoeffs()[patchi] = patchFlux*psf.valueInternalCoeffs(pw);
        fvm.boundaryCoeffs()[patchi] = -patchFlux*psf.valueBoundaryCoeffs(pw);
    }

    // Gradient as linearUpwindV::correction gets it (after the matrix)
    tmp<volTensorField> tgradVf =
        fv::gradScheme<vector>::New(mesh, mesh.gradScheme(gradName))()
       .grad(vf, gradName);
    const volTensorField& gradVf = tgradVf();

    const scalarField& phi = faceFlux.primitiveField();
    const scalarField& w = mesh.weights().primitiveField();
    const vectorField& C = mesh.C().primitiveField();
    const vectorField& Cf = mesh.Cf().primitiveField();
    const vectorField& v = vf.primitiveField();
    const tensorField& gv = gradVf.primitiveField();

    vectorField acc(vf.size(), Zero);

    forAll(own, facei)
    {
        // gaussConvectionScheme: upwind coefficients
        const scalar wf = pos0(phi[facei]);
        const scalar lo = (-wf)*phi[facei];
        const scalar up = lo + phi[facei];
        lower[facei] = lo;
        upper[facei] = up;
        diag[own[facei]] -= lo;
        diag[nei[facei]] -= up;

        // linearUpwindV::correction
        vector maxCorr;
        vector sfCorr;

        if (phi[facei] > 0.0)
        {
            maxCorr = (1.0 - w[facei])*(v[nei[facei]] - v[own[facei]]);
            sfCorr = (Cf[facei] - C[own[facei]]) & gv[own[facei]];
        }
        else
        {
            maxCorr = w[facei]*(v[own[facei]] - v[nei[facei]]);
            sfCorr = (Cf[facei] - C[nei[facei]]) & gv[nei[facei]];
        }

        const scalar sfCorrs = magSqr(sfCorr);
        const scalar maxCorrs = sfCorr & maxCorr;

        if (sfCorrs > 0)
        {
            if (maxCorrs < 0)
            {
                sfCorr = Zero;
            }
            else if (sfCorrs > maxCorrs)
            {
                sfCorr *= maxCorrs/(sfCorrs + VSMALL);
            }
        }

        // faceFlux*correction, surfaceIntegrate
        const vector flux = phi[facei]*sfCorr;
        acc[own[facei]] += flux;
        acc[nei[facei]] -= flux;
    }

    forAll(mesh.boundary(), patchi)
    {
        const labelUList& pOwner = mesh.boundary()[patchi].faceCells();
        const scalarField& pFaceFlux = faceFlux.boundaryField()[patchi];

        // Zero correction except on coupled patches (as linearUpwindV)
        vectorField pSfCorr(pOwner.size(), Zero);

        if (vf.boundaryField()[patchi].coupled())
        {
            const vectorField& pCf = mesh.Cf().boundaryField()[patchi];
            const scalarField& pW = mesh.weights().boundaryField()[patchi];

            const tensorField pGradVfNei
            (
                gradVf.boundaryField()[patchi].patchNeighbourField()
            );
            const vectorField pVfNei
            (
                vf.boundaryField()[patchi].patchNeighbourField()
            );
            const vectorField pd(mesh.Cf().boundaryField()[patchi].patch().delta());

            forAll(pOwner, facei)
            {
                const label o = pOwner[facei];

                vector maxCorr;

                if (pFaceFlux[facei] > 0)
                {
                    pSfCorr[facei] = (pCf[facei] - C[o]) & gv[o];

                    maxCorr = (1.0 - pW[facei])*(pVfNei[facei] - v[o]);
                }
                else
                {
                    pSfCorr[facei] =
                        (pCf[facei] - pd[facei] - C[o]) & pGradVfNei[facei];

                    maxCorr = pW[facei]*(v[o] - pVfNei[facei]);
                }

                const scalar pSfCorrs = magSqr(pSfCorr[facei]);
                const scalar maxCorrs = pSfCorr[facei] & maxCorr;

                if (pSfCorrs > 0)
                {
                    if (maxCorrs < 0)
                    {
                        pSfCorr[facei] = Zero;
                    }
                    else if (pSfCorrs > maxCorrs)
                    {
                        pSfCorr[facei] *= maxCorrs/(pSfCorrs + VSMALL);
                    }
                }
            }
        }

        forAll(pOwner, facei)
        {
            acc[pOwner[facei]] += pFaceFlux[facei]*pSfCorr[facei];
        }
    }

    acc /= mesh.Vsc();
    fvm.source() -= mesh.V()*acc;

    return tfvm;
}

} // End namespace Foam


// ************************************************************************* //
