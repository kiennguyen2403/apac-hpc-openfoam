/*---------------------------------------------------------------------------*\
    simpleFoamFused: steps of the simpleFoam momentum and pressure equations
    computed in fewer passes over memory, with the same values as the stock
    code (bit-identical).

    1. solve(UEqn == -fvc::grad(p)) copies the whole momentum matrix to add
       V*(-grad(p)) to its source. Here the source with that term is built
       separately and swapped into UEqn for the solve, then swapped back
       (the solve leaves the matrix unchanged).
    2. phiHbyA = flux(HbyA) [+ makeRelative, adjustPhi: no-ops here] and the
       SIMPLEC term interpolate(rAtU - rAU)*snGrad(p)*magSf: one face loop
       instead of the flux, interpolation, snGrad, correction and three
       product/sum face fields.
    3. phi = phiHbyA - pEqn.flux(): one face loop instead of the flux field
       (faceH, + face flux correction) and the difference field.
    4. HbyA -= (rAU - rAtU)*grad(p) and U = HbyA - rAtU*grad(p): one cell
       loop each instead of two temporary fields.

    Each value is the same expression as in the stock code: the linear
    interpolation lambda*(a - b) + b and faceH's Upper*psi_u - Lower*psi_l
    are written as there (same contraction); every other operation that
    the stock code does as a separate field operation is a separate
    statement here (not contracted). Patch values use the stock functions.

    SIMPLEFOAM_FUSED_OFF=1 runs the stock code; SIMPLEFOAM_FUSED_REPORT=1
    prints (stderr, master, once) which parts ran.
\*---------------------------------------------------------------------------*/

#ifndef simpleFoamFusedSteps_H
#define simpleFoamFusedSteps_H

#include "fvCFD.H"
#include "HashSet.H"

#include <cstdlib>
#include <iostream>

namespace Foam
{
namespace simpleFoamFused
{

inline bool off()
{
    return std::getenv("SIMPLEFOAM_FUSED_OFF");
}


inline void report(const char* what)
{
    static wordHashSet done;

    if
    (
        std::getenv("SIMPLEFOAM_FUSED_REPORT")
     && Pstream::master()
     && done.insert(what)
    )
    {
        std::cerr << "simpleFoamFused: " << what << std::endl;
    }
}


inline bool schemeIs(const ITstream& is, const char* name)
{
    return is.size() == 1 && is[0].isWord() && is[0].wordToken() == name;
}


// The pressure-equation steps 2-4 reproduce the stock code only for
// SIMPLEC (consistent), no MRF zones, adjustPhi a no-op (p has a fixed
// value), linear interpolation and the corrected snGrad
inline bool pEqnOk
(
    const fvMesh& mesh,
    const volScalarField& p,
    const label nMRF,
    const bool consistent,
    const volVectorField& HbyA,
    const volScalarField& d
)
{
    return
        !off()
     && consistent
     && nMRF == 0
     && !p.needReference()
     && schemeIs(mesh.interpolationScheme("flux(" + HbyA.name() + ')'), "linear")
     && schemeIs(mesh.interpolationScheme("interpolate(" + d.name() + ')'), "linear")
     && schemeIs(mesh.snGradScheme("snGrad(" + p.name() + ')'), "corrected");
}


// 1. UEqn source + V*(-grad(p)), as operator==(fvMatrix, tmp<field>)
inline vectorField momentumSource
(
    const fvVectorMatrix& UEqn,
    const volVectorField& gradp
)
{
    vectorField b(UEqn.source());

    const scalarField& V = gradp.mesh().V();
    const vectorField& g = gradp.primitiveField();

    forAll(b, i)
    {
        const vector ng = -g[i];
        const vector Vng = V[i]*ng;
        b[i] += Vng;
    }

    return b;
}


// 2. flux(HbyA) + interpolate(d)*snGrad(p)*magSf, d = rAtU - rAU
inline tmp<surfaceScalarField> phiHbyA
(
    const volVectorField& HbyA,
    const volScalarField& d,
    const volScalarField& p
)
{
    const fvMesh& mesh = p.mesh();
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();

    const surfaceScalarField& lambdas = mesh.surfaceInterpolation::weights();
    const surfaceVectorField& Sf = mesh.Sf();
    const surfaceScalarField& magSf = mesh.magSf();
    const surfaceScalarField& deltaCoeffs = mesh.nonOrthDeltaCoeffs();
    const surfaceVectorField& ncv = mesh.nonOrthCorrectionVectors();

    // grad(p) as correctedSnGrad::fullGradCorrection gets it
    tmp<volVectorField> tgradp =
        fv::gradScheme<scalar>::New
        (
            mesh,
            mesh.gradScheme("grad(" + p.name() + ')')
        )().grad(p, "grad(" + p.name() + ')');
    const volVectorField& gradp = tgradp();

    tmp<surfaceScalarField> tphi
    (
        new surfaceScalarField
        (
            IOobject
            (
                "interpolate(" + HbyA.name() + ')',
                HbyA.instance(),
                HbyA.db()
            ),
            mesh,
            Sf.dimensions()*HbyA.dimensions()
        )
    );
    surfaceScalarField& phi = tphi.ref();
    phi.setOriented();

    const scalarField& lambda = lambdas.primitiveField();
    const vectorField& Sfi = Sf.primitiveField();
    const scalarField& mSf = magSf.primitiveField();
    const scalarField& dC = deltaCoeffs.primitiveField();
    const vectorField& ncvi = ncv.primitiveField();
    const vectorField& H = HbyA.primitiveField();
    const scalarField& di = d.primitiveField();
    const scalarField& pv = p.primitiveField();
    const vectorField& G = gradp.primitiveField();
    scalarField& phii = phi.primitiveFieldRef();

    forAll(own, facei)
    {
        const label o = own[facei];
        const label n = nei[facei];

        // fvc::flux(HbyA)
        const scalar flux =
            Sfi[facei] & (lambda[facei]*(H[o] - H[n]) + H[n]);

        // fvc::interpolate(rAtU - rAU)
        const scalar dface = lambda[facei]*(di[o] - di[n]) + di[n];

        // fvc::snGrad(p): deltaCoeffs*(p_N - p_P), then += correction
        scalar sn = dC[facei]*(pv[n] - pv[o]);
        const scalar corr =
            ncvi[facei] & (lambda[facei]*(G[o] - G[n]) + G[n]);
        sn += corr;

        // interpolate*snGrad, *magSf, phiHbyA +=
        const scalar t2 = dface*sn;
        const scalar t3 = t2*mSf[facei];
        phii[facei] = flux + t3;
    }

    forAll(mesh.boundary(), patchi)
    {
        const fvsPatchScalarField& pLambda = lambdas.boundaryField()[patchi];
        const fvsPatchVectorField& pSf = Sf.boundaryField()[patchi];

        const fvPatchVectorField& pH = HbyA.boundaryField()[patchi];
        scalarField pflux;
        if (pH.coupled())
        {
            pflux =
                pSf
              & lerp
                (
                    pH.patchNeighbourField(),
                    pH.patchInternalField(),
                    pLambda
                );
        }
        else
        {
            pflux = pSf & pH;
        }

        const fvPatchScalarField& pd = d.boundaryField()[patchi];
        scalarField pdface;
        if (pd.coupled())
        {
            pdface =
                lerp(pd.patchNeighbourField(), pd.patchInternalField(), pLambda);
        }
        else
        {
            pdface = pd;
        }

        const fvPatchScalarField& pp = p.boundaryField()[patchi];
        scalarField psn;
        if (pp.coupled())
        {
            psn = pp.snGrad(deltaCoeffs.boundaryField()[patchi]);
        }
        else
        {
            psn = pp.snGrad();
        }

        const fvPatchVectorField& pG = gradp.boundaryField()[patchi];
        const fvsPatchVectorField& pncv = ncv.boundaryField()[patchi];
        scalarField pcorr;
        if (pG.coupled())
        {
            pcorr =
                pncv
              & lerp(pG.patchNeighbourField(), pG.patchInternalField(), pLambda);
        }
        else
        {
            pcorr = pncv & pG;
        }
        psn += pcorr;

        scalarField& pphi = phi.boundaryFieldRef()[patchi];
        pphi = pflux;
        pphi += (pdface*psn)*magSf.boundaryField()[patchi];
    }

    return tphi;
}


// 3. phiHbyA - pEqn.flux() (flux: faceH + face flux correction)
inline tmp<surfaceScalarField> phiCorrected
(
    const surfaceScalarField& phiHbyA,
    fvScalarMatrix& pEqn,
    const volScalarField& p
)
{
    const fvMesh& mesh = p.mesh();
    const lduMatrix& m = pEqn;
    const scalarField& Lower = m.lower();
    const scalarField& Upper = m.upper();
    const labelUList& l = m.lduAddr().lowerAddr();
    const labelUList& u = m.lduAddr().upperAddr();
    const scalarField& psi = p.primitiveField();

    tmp<surfaceScalarField> tres
    (
        new surfaceScalarField
        (
            IOobject
            (
                "phiCorrected",
                p.instance(),
                p.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh,
            phiHbyA.dimensions()
        )
    );
    surfaceScalarField& res = tres.ref();
    res.setOriented();

    const scalarField& ph = phiHbyA.primitiveField();
    scalarField& r = res.primitiveFieldRef();

    const auto& ffcPtr = pEqn.faceFluxCorrectionPtr();

    if (ffcPtr)
    {
        const scalarField& ffc = ffcPtr->primitiveField();
        for (label face=0; face<l.size(); face++)
        {
            const scalar fH =
                Upper[face]*psi[u[face]]
              - Lower[face]*psi[l[face]];
            const scalar fl = fH + ffc[face];
            r[face] = ph[face] - fl;
        }
    }
    else
    {
        for (label face=0; face<l.size(); face++)
        {
            const scalar fH =
                Upper[face]*psi[u[face]]
              - Lower[face]*psi[l[face]];
            r[face] = ph[face] - fH;
        }
    }

    forAll(mesh.boundary(), patchi)
    {
        const fvPatchScalarField& pp = p.boundaryField()[patchi];

        const scalarField ic
        (
            cmptMultiply
            (
                pEqn.internalCoeffs()[patchi],
                pp.patchInternalField()
            )
        );

        scalarField nc(pEqn.boundaryCoeffs()[patchi]);
        if (pp.coupled())
        {
            nc = cmptMultiply
            (
                pEqn.boundaryCoeffs()[patchi],
                pp.patchNeighbourField()
            );
        }

        scalarField fb(ic - nc);
        if (ffcPtr)
        {
            fb += ffcPtr->boundaryField()[patchi];
        }

        res.boundaryFieldRef()[patchi] = phiHbyA.boundaryField()[patchi] - fb;
    }

    return tres;
}


// 4a. HbyA -= (rAU - rAtU)*grad(p)
inline void HbyAUpdate
(
    volVectorField& HbyA,
    const volScalarField& rAU,
    const volScalarField& rAtU,
    const volVectorField& gradp
)
{
    vectorField& h = HbyA.primitiveFieldRef();
    const scalarField& a = rAU.primitiveField();
    const scalarField& b = rAtU.primitiveField();
    const vectorField& g = gradp.primitiveField();

    forAll(h, i)
    {
        const scalar t = a[i] - b[i];
        const vector v = t*g[i];
        h[i] -= v;
    }

    auto& hbf = HbyA.boundaryFieldRef();
    forAll(hbf, patchi)
    {
        hbf[patchi] -=
        (
            (rAU.boundaryField()[patchi] - rAtU.boundaryField()[patchi])
           *gradp.boundaryField()[patchi]
        );
    }
}


// 4b. HbyA - rAtU*grad(p), as a new field (assigned to U by the caller)
inline tmp<volVectorField> correctedVelocity
(
    const volVectorField& HbyA,
    const volScalarField& rAtU,
    const volVectorField& gradp
)
{
    const fvMesh& mesh = HbyA.mesh();

    tmp<volVectorField> tres
    (
        new volVectorField
        (
            IOobject
            (
                "UCorrected",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh,
            dimensioned<vector>(HbyA.dimensions(), Zero)
        )
    );
    volVectorField& res = tres.ref();

    vectorField& r = res.primitiveFieldRef();
    const vectorField& h = HbyA.primitiveField();
    const scalarField& b = rAtU.primitiveField();
    const vectorField& g = gradp.primitiveField();

    forAll(r, i)
    {
        const vector v = b[i]*g[i];
        r[i] = h[i] - v;
    }

    forAll(res.boundaryField(), patchi)
    {
        vectorField& rp = res.boundaryFieldRef()[patchi];
        rp =
            HbyA.boundaryField()[patchi]
          - rAtU.boundaryField()[patchi]*gradp.boundaryField()[patchi];
    }

    return tres;
}

} // End namespace simpleFoamFused
} // End namespace Foam

#endif
