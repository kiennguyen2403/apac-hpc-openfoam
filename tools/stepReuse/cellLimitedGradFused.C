/*---------------------------------------------------------------------------*\
    Replacement for fv::cellLimitedGrad<vector, minmod>::calcGrad (grad(U)
    cellLimited Gauss linear 1): the in-place version (solution cache update)
    and the returning version (which calls it).

    Stock: the Gauss gradient (interpolate, face loop, patches, /V, boundary
    conditions); then copies of U as max and min, a face loop and patch loop
    taking the neighbour max/min, max -= U and min -= U (two passes), a
    limiter field set to one, a face loop and patch loop applying the minmod
    limiter on both sides of each face, and the gradient scaled by it.

    Here the neighbour max/min is taken in the Gauss gradient's face loop
    (which already reads both cell values), max, min and the limiter are
    initialised in one loop, and the deltas max - U, min - U are formed where
    the limiter loop uses them (the same subtraction). The Gauss sums keep
    their order; max, min and the limiter are min/max selections, exact in
    any order. Bit-identical.

    Only for cellLimited Gauss linear with coefficient 1; otherwise, or with
    STEP_REUSE_OFF=climit, the stock function runs.
\*---------------------------------------------------------------------------*/

#include "cellLimitedGrad.H"
#include "minmodGradientLimiter.H"
#include "gaussGrad.H"
#include "stepReuse.H"

namespace Foam
{

typedef fv::cellLimitedGrad<vector, fv::gradientLimiters::minmod>
    stepReuseCellLimitedVector;

template<>
void stepReuseCellLimitedVector::calcGrad
(
    volTensorField& g,
    const volVectorField& vsf
) const
{
    const fvMesh& mesh = vsf.mesh();

    bool fused = !stepReuse::off("climit") && k_ >= 1.0;

    if (fused)
    {
        const ITstream& is = mesh.gradScheme(g.name());
        fused =
            is.size() == 4
         && is[0].isWord() && is[0].wordToken() == "cellLimited"
         && is[1].isWord() && is[1].wordToken() == "Gauss"
         && is[2].isWord() && is[2].wordToken() == "linear";
    }

    if (!fused)
    {
        typedef void (*fn)
        (
            const stepReuseCellLimitedVector*,
            volTensorField&,
            const volVectorField&
        );
        static const fn stock = reinterpret_cast<fn>
        (
            stepReuse::next
            (
                "_ZNK4Foam2fv15cellLimitedGradINS_6VectorIdEENS0_16gradientLimiters"
                "6minmodEE8calcGradERNS_14GeometricFieldINS_6TensorIdEENS_"
                "12fvPatchFieldENS_7volMeshEEERKNS7_IS3_SA_SB_EE"
            )
        );
        stock(this, g, vsf);
        return;
    }

    stepReuse::report("climit", "cellLimitedGrad: limiter passes fused");

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const surfaceScalarField& lambdas = mesh.surfaceInterpolation::weights();
    const scalarField& lambda = lambdas.primitiveField();
    const vectorField& Sf = mesh.Sf().primitiveField();
    const vectorField& C = mesh.C().primitiveField();
    const surfaceVectorField& Cf = mesh.Cf();
    const vectorField& Cfi = Cf.primitiveField();

    const vectorField& v = vsf.primitiveField();
    const label nCells = v.size();

    // max, min (start: U) and the limiter (start: one)
    vectorField maxVsf(nCells);
    vectorField minVsf(nCells);
    vectorField limiter(nCells);
    for (label c = 0; c < nCells; ++c)
    {
        maxVsf[c] = v[c];
        minVsf[c] = v[c];
        limiter[c] = pTraits<vector>::one;
    }

    // Gauss gradient (gaussGrad::calcGrad in place) + neighbour max/min
    g = dimensioned<tensor>(vsf.dimensions()/dimLength, Zero);
    tensorField& ig = g.primitiveFieldRef();

    forAll(owner, facei)
    {
        const label own = owner[facei];
        const label nei = neighbour[facei];

        const vector ssf = lambda[facei]*(v[own] - v[nei]) + v[nei];
        const tensor Sfssf = Sf[facei]*ssf;

        ig[own] += Sfssf;
        ig[nei] -= Sfssf;

        maxVsf[own] = max(maxVsf[own], v[nei]);
        minVsf[own] = min(minVsf[own], v[nei]);

        maxVsf[nei] = max(maxVsf[nei], v[own]);
        minVsf[nei] = min(minVsf[nei], v[own]);
    }

    forAll(mesh.boundary(), patchi)
    {
        const labelUList& pFaceCells = mesh.boundary()[patchi].faceCells();
        const vectorField& pSf = mesh.Sf().boundaryField()[patchi];
        const fvPatchVectorField& pvf = vsf.boundaryField()[patchi];

        if (pvf.coupled())
        {
            const vectorField pssf
            (
                lerp
                (
                    pvf.patchNeighbourField(),
                    pvf.patchInternalField(),
                    lambdas.boundaryField()[patchi]
                )
            );

            forAll(mesh.boundary()[patchi], facei)
            {
                ig[pFaceCells[facei]] += pSf[facei]*pssf[facei];
            }
        }
        else
        {
            forAll(mesh.boundary()[patchi], facei)
            {
                ig[pFaceCells[facei]] += pSf[facei]*pvf[facei];
            }
        }
    }

    ig /= mesh.V();

    g.correctBoundaryConditions();
    gaussGrad<vector>::correctBoundaryConditions(vsf, g);

    // Neighbour max/min across the patches
    const auto& bsf = vsf.boundaryField();

    forAll(bsf, patchi)
    {
        const fvPatchVectorField& psf = bsf[patchi];
        const labelUList& pOwner = mesh.boundary()[patchi].faceCells();

        if (psf.coupled())
        {
            const vectorField psfNei(psf.patchNeighbourField());

            forAll(pOwner, pFacei)
            {
                const label own = pOwner[pFacei];
                const vector& vsfNei = psfNei[pFacei];

                maxVsf[own] = max(maxVsf[own], vsfNei);
                minVsf[own] = min(minVsf[own], vsfNei);
            }
        }
        else
        {
            forAll(pOwner, pFacei)
            {
                const label own = pOwner[pFacei];
                const vector& vsfNei = psf[pFacei];

                maxVsf[own] = max(maxVsf[own], vsfNei);
                minVsf[own] = min(minVsf[own], vsfNei);
            }
        }
    }

    // minmod limiter on both sides of each face (deltas max - U, min - U)
    forAll(owner, facei)
    {
        const label own = owner[facei];
        const label nei = neighbour[facei];

        limitFace
        (
            limiter[own],
            maxVsf[own] - v[own],
            minVsf[own] - v[own],
            (Cfi[facei] - C[own]) & ig[own]
        );

        limitFace
        (
            limiter[nei],
            maxVsf[nei] - v[nei],
            minVsf[nei] - v[nei],
            (Cfi[facei] - C[nei]) & ig[nei]
        );
    }

    forAll(bsf, patchi)
    {
        const labelUList& pOwner = mesh.boundary()[patchi].faceCells();
        const vectorField& pCf = Cf.boundaryField()[patchi];

        forAll(pOwner, pFacei)
        {
            const label own = pOwner[pFacei];

            limitFace
            (
                limiter[own],
                maxVsf[own] - v[own],
                minVsf[own] - v[own],
                ((pCf[pFacei] - C[own]) & ig[own])
            );
        }
    }

    if (fv::debug)
    {
        auto limits = gMinMax(limiter);
        auto avg = gAverage(limiter);

        Info<< "gradient limiter for: " << vsf.name()
            << " min = " << limits.min()
            << " max = " << limits.max()
            << " average: " << avg << endl;
    }

    limitGradient(limiter, ig);
    g.correctBoundaryConditions();
    gaussGrad<vector>::correctBoundaryConditions(vsf, g);
}


template<>
tmp<volTensorField> stepReuseCellLimitedVector::calcGrad
(
    const volVectorField& vsf,
    const word& name
) const
{
    const fvMesh& mesh = vsf.mesh();

    tmp<volTensorField> tGrad
    (
        new volTensorField
        (
            IOobject
            (
                name,
                vsf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            vsf.dimensions()/dimLength,
            fvPatchFieldBase::extrapolatedCalculatedType()
        )
    );

    calcGrad(tGrad.ref(), vsf);

    return tGrad;
}

} // End namespace Foam


// ************************************************************************* //
