/*---------------------------------------------------------------------------*\
    Replacement for kOmegaSSTBase<incompressible RAS eddyViscosity>::correct()
    (explicit specialisation; the library's weak instantiation and its
    this-adjusting thunk both reach it through the PLT under LD_PRELOAD).

    The v2512 code computes some results more than once per call:

    1. fvc::surfaceIntegrate(phi) three times: for divU, and in the 'bounded'
       correction of both convection terms (bounded Gauss upwind for omega
       and k). Computed once here.

    2. The upwind convection coefficients of omega and k: the weights are
       pos0(phi), independent of the field, so lower, upper and the diagonal
       (after the 'bounded' correction) are identical. Built once for omega,
       kept and moved into the k matrix; only the boundary coefficients,
       which depend on the field's patch types, are built per field.

    (grad(k) and grad(omega), computed for CDkOmega and again by the
    corrected laplacian terms, are reused by the stock solution cache once
    they are listed in fvSolution: cache { grad(k); grad(omega); })

    4. (dev, STEP_REUSE_OFF=kfuse to switch off) The per-cell maths: S2,
       GbyNu0, CDkOmega, F1, F23 (twice: here and in correctNut), GbyNu,
       gamma, beta, DkEff, DomegaEff and nut are each computed in one loop
       over the cells and one per patch, instead of one field operation (and
       one temporary field) per arithmetic operation. Each value goes through
       the same scalar operations in the same order as in the field
       operations, with FMA contraction off in these loops (the stock field
       operations are one operation per loop, so nothing is contracted
       there); the tensor functions (symm, magSqr, devTwoSymm, &&, &) are
       the same inline functions. Only for the plain kOmegaSST model (type
       kOmegaSST) with F3 off.

    Every result is produced by the same operations on the same values as in
    the stock code: bit-identical. Anything else (other schemes, moving mesh,
    STEP_REUSE_OFF=komega) runs the stock correct().
\*---------------------------------------------------------------------------*/

#include "IncompressibleTurbulenceModel.H"
#include "incompressible/transportModel/transportModel.H"
#include "RASModel.H"
#include "eddyViscosity.H"
#include "kOmegaSSTBase.H"
#include "fvOptions.H"
#include "bound.H"
#include "fvcSurfaceIntegrate.H"
#include "stepReuse.H"

namespace Foam
{

typedef eddyViscosity<RASModel<IncompressibleTurbulenceModel<transportModel>>>
    stepReuseEddyViscosityRAS;

namespace stepReuseDetail
{

// The div scheme is exactly "bounded Gauss upwind"
static bool boundedGaussUpwind(const ITstream& is)
{
    return
        is.size() == 3
     && is[0].isWord() && is[0].wordToken() == "bounded"
     && is[1].isWord() && is[1].wordToken() == "Gauss"
     && is[2].isWord() && is[2].wordToken() == "upwind";
}


// gaussConvectionScheme<scalar>::fvmDiv with upwind weights, optionally
// taking lower/upper/diag from an earlier field (moved, not copied)
static tmp<fvScalarMatrix> upwindDiv
(
    const surfaceScalarField& faceFlux,
    const surfaceScalarField& weights,
    const volScalarField& vf,
    scalarField* lowerPtr = nullptr,
    scalarField* upperPtr = nullptr,
    scalarField* diagPtr = nullptr
)
{
    tmp<fvScalarMatrix> tfvm
    (
        new fvScalarMatrix
        (
            vf,
            faceFlux.dimensions()*vf.dimensions()
        )
    );
    fvScalarMatrix& fvm = tfvm.ref();

    if (lowerPtr)
    {
        fvm.lower().transfer(*lowerPtr);
        fvm.upper().transfer(*upperPtr);
        fvm.diag().transfer(*diagPtr);
    }
    else
    {
        fvm.lower() = -weights.primitiveField()*faceFlux.primitiveField();
        fvm.upper() = fvm.lower() + faceFlux.primitiveField();
        fvm.negSumDiag();
    }

    forAll(vf.boundaryField(), patchi)
    {
        const fvPatchField<scalar>& psf = vf.boundaryField()[patchi];
        const fvsPatchScalarField& patchFlux = faceFlux.boundaryField()[patchi];
        const fvsPatchScalarField& pw = weights.boundaryField()[patchi];

        fvm.internalCoeffs()[patchi] = patchFlux*psf.valueInternalCoeffs(pw);
        fvm.boundaryCoeffs()[patchi] = -patchFlux*psf.valueBoundaryCoeffs(pw);
    }

    return tfvm;
}



// ---------------------------------------------------------------------------
// 4. Fused per-cell maths (see the header). Constants are the values of the
//    dimensioned products the stock expressions form.

struct sstConstants
{
    scalar cA1;     // scalar(1)/betaStar
    scalar cA2;     // scalar(2)/betaStar
    scalar c2aO2;   // 2*alphaOmega2
    scalar c4aO2;   // 4*alphaOmega2
    scalar cG;      // (c1/a1)*betaStar
    scalar a1;
    scalar b1;
};


// New unregistered field; calculated patches (constraint types on constraint
// patches), as the stock field operations make them
static tmp<volScalarField> newScalarField
(
    const word& name,
    const fvMesh& mesh,
    const dimensionSet& dims
)
{
    return tmp<volScalarField>::New
    (
        IOobject
        (
            name,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh,
        dimensionedScalar(dims, Zero)
    );
}

static tmp<volScalarField::Internal> newInternalField
(
    const word& name,
    const fvMesh& mesh,
    const dimensionSet& dims
)
{
    return tmp<volScalarField::Internal>::New
    (
        IOobject
        (
            name,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh,
        dimensionedScalar(dims, Zero)
    );
}

// 2*magSqr(symm(gradU))
static inline scalar S2Of(const tensor& g)
{
    return 2*magSqr(symm(g));
}

// (2*alphaOmega2)*(gradk & gradomega)/omega
static inline scalar CDkOmegaOf
(
    const vector& gk, const vector& go, const scalar o, const sstConstants& c
)
{
    #pragma clang fp contract(off)
    return c.c2aO2*(gk & go)/o;
}

// F1: tanh(pow4(min(min(max(A, B), C), 10)))
static inline scalar F1Of
(
    const scalar k, const scalar o, const scalar y, const scalar nu,
    const scalar cd, const sstConstants& c
)
{
    #pragma clang fp contract(off)
    const scalar cdp = max(cd, scalar(1.0e-10));
    const scalar A = c.cA1*sqrt(k)/(o*y);
    const scalar B = scalar(500)*nu/(sqr(y)*o);
    const scalar C = c.c4aO2*k/(cdp*sqr(y));
    const scalar arg1 = min(min(max(A, B), C), scalar(10));
    return tanh(pow4(arg1));
}

// F2: tanh(sqr(min(max(A2, B), 100)))
static inline scalar F2Of
(
    const scalar k, const scalar o, const scalar y, const scalar nu,
    const sstConstants& c
)
{
    #pragma clang fp contract(off)
    const scalar A = c.cA2*sqrt(k)/(o*y);
    const scalar B = scalar(500)*nu/(sqr(y)*o);
    const scalar arg2 = min(max(A, B), scalar(100));
    return tanh(sqr(arg2));
}

// F1*(psi1 - psi2) + psi2
static inline scalar blendOf(const scalar f1, const scalar d, const scalar p2)
{
    #pragma clang fp contract(off)
    const scalar t = f1*d;
    return t + p2;
}

// blend(F1, psi1, psi2)*nut + nu
static inline scalar DEffOf
(
    const scalar f1, const scalar d, const scalar p2,
    const scalar nut, const scalar nu
)
{
    #pragma clang fp contract(off)
    const scalar a = blendOf(f1, d, p2);
    const scalar t = a*nut;
    return t + nu;
}

// a1*k/max(a1*omega, b1*F23*sqrt(S2))
static inline scalar nutOf
(
    const scalar k, const scalar o, const scalar f23, const scalar s2,
    const sstConstants& c
)
{
    #pragma clang fp contract(off)
    return c.a1*k/max(c.a1*o, c.b1*f23*sqrt(s2));
}

// min(GbyNu0, (c1/a1)*betaStar*omega*max(a1*omega, b1*F2*sqrt(S2)))
static inline scalar GbyNuOf
(
    const scalar g0, const scalar o, const scalar f2, const scalar s2,
    const sstConstants& c
)
{
    #pragma clang fp contract(off)
    return min(g0, c.cG*o*max(c.a1*o, c.b1*f2*sqrt(s2)));
}


static tmp<volScalarField> fusedS2(const volTensorField& gradU)
{
    tmp<volScalarField> tres =
        newScalarField("S2", gradU.mesh(), sqr(gradU.dimensions()));
    volScalarField& res = tres.ref();

    scalarField& r = res.primitiveFieldRef();
    const tensorField& g = gradU.primitiveField();
    forAll(r, i) { r[i] = S2Of(g[i]); }

    forAll(res.boundaryField(), pi)
    {
        scalarField& rp = res.boundaryFieldRef()[pi];
        const tensorField& gp = gradU.boundaryField()[pi];
        forAll(rp, i) { rp[i] = S2Of(gp[i]); }
    }
    return tres;
}


static tmp<volScalarField::Internal> fusedGbyNu0
(
    const word& name,
    const volTensorField& gradU
)
{
    tmp<volScalarField::Internal> tres =
        newInternalField(name, gradU.mesh(), sqr(gradU.dimensions()));
    scalarField& r = tres.ref().field();
    const tensorField& g = gradU.primitiveField();
    forAll(r, i) { r[i] = g[i] && devTwoSymm(g[i]); }
    return tres;
}


static tmp<volScalarField> fusedCDkOmega
(
    const volVectorField& gk,
    const volVectorField& go,
    const volScalarField& omega,
    const sstConstants& c
)
{
    tmp<volScalarField> tres = newScalarField
    (
        "CDkOmega",
        omega.mesh(),
        gk.dimensions()*go.dimensions()/omega.dimensions()
    );
    volScalarField& res = tres.ref();

    scalarField& r = res.primitiveFieldRef();
    const vectorField& gki = gk.primitiveField();
    const vectorField& goi = go.primitiveField();
    const scalarField& oi = omega.primitiveField();
    forAll(r, i) { r[i] = CDkOmegaOf(gki[i], goi[i], oi[i], c); }

    forAll(res.boundaryField(), pi)
    {
        scalarField& rp = res.boundaryFieldRef()[pi];
        const vectorField& gkp = gk.boundaryField()[pi];
        const vectorField& gop = go.boundaryField()[pi];
        const scalarField& op = omega.boundaryField()[pi];
        forAll(rp, i) { rp[i] = CDkOmegaOf(gkp[i], gop[i], op[i], c); }
    }
    return tres;
}


static tmp<volScalarField> fusedF1
(
    const volScalarField& k, const volScalarField& omega,
    const volScalarField& y, const volScalarField& nu,
    const volScalarField& cd, const sstConstants& c
)
{
    tmp<volScalarField> tres = newScalarField("F1", k.mesh(), dimless);
    volScalarField& res = tres.ref();

    scalarField& r = res.primitiveFieldRef();
    const scalarField& ki = k.primitiveField();
    const scalarField& oi = omega.primitiveField();
    const scalarField& yi = y.primitiveField();
    const scalarField& nui = nu.primitiveField();
    const scalarField& cdi = cd.primitiveField();
    forAll(r, i) { r[i] = F1Of(ki[i], oi[i], yi[i], nui[i], cdi[i], c); }

    forAll(res.boundaryField(), pi)
    {
        scalarField& rp = res.boundaryFieldRef()[pi];
        const scalarField& kp = k.boundaryField()[pi];
        const scalarField& op = omega.boundaryField()[pi];
        const scalarField& yp = y.boundaryField()[pi];
        const scalarField& nup = nu.boundaryField()[pi];
        const scalarField& cdp = cd.boundaryField()[pi];
        forAll(rp, i) { rp[i] = F1Of(kp[i], op[i], yp[i], nup[i], cdp[i], c); }
    }
    return tres;
}


static tmp<volScalarField> fusedF2
(
    const volScalarField& k, const volScalarField& omega,
    const volScalarField& y, const volScalarField& nu,
    const sstConstants& c
)
{
    tmp<volScalarField> tres = newScalarField("F23", k.mesh(), dimless);
    volScalarField& res = tres.ref();

    scalarField& r = res.primitiveFieldRef();
    const scalarField& ki = k.primitiveField();
    const scalarField& oi = omega.primitiveField();
    const scalarField& yi = y.primitiveField();
    const scalarField& nui = nu.primitiveField();
    forAll(r, i) { r[i] = F2Of(ki[i], oi[i], yi[i], nui[i], c); }

    forAll(res.boundaryField(), pi)
    {
        scalarField& rp = res.boundaryFieldRef()[pi];
        const scalarField& kp = k.boundaryField()[pi];
        const scalarField& op = omega.boundaryField()[pi];
        const scalarField& yp = y.boundaryField()[pi];
        const scalarField& nup = nu.boundaryField()[pi];
        forAll(rp, i) { rp[i] = F2Of(kp[i], op[i], yp[i], nup[i], c); }
    }
    return tres;
}


// blend(F1, psi1, psi2), cell values (gamma, beta)
static tmp<volScalarField::Internal> fusedBlend
(
    const word& name,
    const volScalarField& F1,
    const dimensionedScalar& psi1,
    const dimensionedScalar& psi2
)
{
    const dimensionedScalar d(psi1 - psi2);
    tmp<volScalarField::Internal> tres =
        newInternalField(name, F1.mesh(), F1.dimensions()*d.dimensions());
    scalarField& r = tres.ref().field();
    const scalarField& f = F1.primitiveField();
    const scalar dv = d.value();
    const scalar p2 = psi2.value();
    forAll(r, i) { r[i] = blendOf(f[i], dv, p2); }
    return tres;
}


// blend(F1, psi1, psi2)*nut + nu (DkEff, DomegaEff)
static tmp<volScalarField> fusedDEff
(
    const word& name,
    const volScalarField& F1,
    const dimensionedScalar& psi1,
    const dimensionedScalar& psi2,
    const volScalarField& nut,
    const volScalarField& nu
)
{
    const dimensionedScalar d(psi1 - psi2);
    tmp<volScalarField> tres = newScalarField(name, F1.mesh(), nut.dimensions());
    volScalarField& res = tres.ref();
    const scalar dv = d.value();
    const scalar p2 = psi2.value();

    scalarField& r = res.primitiveFieldRef();
    const scalarField& f = F1.primitiveField();
    const scalarField& nt = nut.primitiveField();
    const scalarField& n = nu.primitiveField();
    forAll(r, i) { r[i] = DEffOf(f[i], dv, p2, nt[i], n[i]); }

    forAll(res.boundaryField(), pi)
    {
        scalarField& rp = res.boundaryFieldRef()[pi];
        const scalarField& fp = F1.boundaryField()[pi];
        const scalarField& ntp = nut.boundaryField()[pi];
        const scalarField& np = nu.boundaryField()[pi];
        forAll(rp, i) { rp[i] = DEffOf(fp[i], dv, p2, ntp[i], np[i]); }
    }
    return tres;
}


static tmp<volScalarField> fusedNut
(
    const volScalarField& k, const volScalarField& omega,
    const volScalarField& F23, const volScalarField& S2,
    const dimensionSet& dims, const sstConstants& c
)
{
    tmp<volScalarField> tres = newScalarField("nut", k.mesh(), dims);
    volScalarField& res = tres.ref();

    scalarField& r = res.primitiveFieldRef();
    const scalarField& ki = k.primitiveField();
    const scalarField& oi = omega.primitiveField();
    const scalarField& fi = F23.primitiveField();
    const scalarField& si = S2.primitiveField();
    forAll(r, i) { r[i] = nutOf(ki[i], oi[i], fi[i], si[i], c); }

    forAll(res.boundaryField(), pi)
    {
        scalarField& rp = res.boundaryFieldRef()[pi];
        const scalarField& kp = k.boundaryField()[pi];
        const scalarField& op = omega.boundaryField()[pi];
        const scalarField& fp = F23.boundaryField()[pi];
        const scalarField& sp = S2.boundaryField()[pi];
        forAll(rp, i) { rp[i] = nutOf(kp[i], op[i], fp[i], sp[i], c); }
    }
    return tres;
}


// GbyNu0 = GbyNu(GbyNu0, F23, S2), in place (cell values)
static void fusedGbyNu
(
    volScalarField::Internal& GbyNu0,
    const volScalarField& omega,
    const volScalarField& F2,
    const volScalarField& S2,
    const sstConstants& c
)
{
    scalarField& g = GbyNu0.field();
    const scalarField& oi = omega.primitiveField();
    const scalarField& fi = F2.primitiveField();
    const scalarField& si = S2.primitiveField();
    forAll(g, i) { g[i] = GbyNuOf(g[i], oi[i], fi[i], si[i], c); }
}

} // End namespace stepReuseDetail


template<>
void kOmegaSSTBase<stepReuseEddyViscosityRAS>::correct()
{
    using namespace stepReuseDetail;
    typedef stepReuseEddyViscosityRAS BasicEddyViscosityModel;

    if (!this->turbulence_)
    {
        return;
    }

    const fvMesh& mesh = this->mesh_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;

    const bool enabled = !stepReuse::off("komega");

    // divU uses phi(), the convection terms alphaRhoPhi: the same field here
    const bool samePhi = (&this->phi()() == &alphaRhoPhi);

    const bool shareConv =
        enabled
     && !mesh.moving()
     && samePhi
     && boundedGaussUpwind
        (
            mesh.divScheme("div(" + alphaRhoPhi.name() + ',' + omega_.name() + ')')
        )
     && boundedGaussUpwind
        (
            mesh.divScheme("div(" + alphaRhoPhi.name() + ',' + k_.name() + ')')
        );

    // 4. Fused per-cell maths: plain kOmegaSST only (no overridden helpers)
    const bool kfuse =
        shareConv
     && !stepReuse::off("kfuse")
     && this->type() == "kOmegaSST"
     && !F3_;

    if (!shareConv)
    {
        stepReuse::report("komega", "kOmegaSST: stock correct()");

        typedef void (*correctFn)(kOmegaSSTBase<BasicEddyViscosityModel>*);
        static const correctFn stock = reinterpret_cast<correctFn>
        (
            stepReuse::next
            (
                "_ZN4Foam13kOmegaSSTBaseINS_13eddyViscosityINS_8RASModelINS_"
                "29IncompressibleTurbulenceModelINS_14transportModelEEEEEEEE"
                "7correctEv"
            )
        );
        stock(this);
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    BasicEddyViscosityModel::correct();

    // 1. Flux divergence for divU and both 'bounded' convection corrections
    //    (stock: fvc::div(fvc::absolute(phi(), U)), the same values for a
    //    static mesh)
    tmp<volScalarField> tdivPhi(fvc::surfaceIntegrate(alphaRhoPhi));
    const volScalarField::Internal& divU = tdivPhi();

    const sstConstants c
    {
        (scalar(1)/betaStar_).value(),
        (scalar(2)/betaStar_).value(),
        (2*alphaOmega2_).value(),
        (4*alphaOmega2_).value(),
        ((c1_/a1_)*betaStar_).value(),
        a1_.value(),
        b1_.value()
    };

    // Laminar viscosity as F1/F2 (mu()/rho = nu()) and DkEff/DomegaEff (nu())
    // use it
    tmp<volScalarField> tnu;
    if (kfuse)
    {
        tnu = this->nu();
    }

    tmp<volTensorField> tgradU = fvc::grad(U);
    const volScalarField S2
    (
        kfuse ? fusedS2(tgradU()) : this->S2(tgradU())
    );
    volScalarField::Internal GbyNu0
    (
        kfuse
      ? fusedGbyNu0(IOobject::scopedName(this->type(), "GbyNu"), tgradU())
      : this->GbyNu0(tgradU(), S2)
    );
    volScalarField::Internal G(this->GName(), nut*GbyNu0);

    // Update omega and G at the wall (see kOmegaSSTBase.C)
    omega_.boundaryFieldRef().updateCoeffs();
    omega_.boundaryFieldRef().template evaluateCoupled<coupledFvPatch>();

    const volScalarField CDkOmega
    (
        kfuse
      ? fusedCDkOmega(fvc::grad(k_)(), fvc::grad(omega_)(), omega_, c)
      : (2*alphaOmega2_)*(fvc::grad(k_) & fvc::grad(omega_))/omega_
    );

    const volScalarField F1
    (
        kfuse
      ? fusedF1(k_, omega_, y_, tnu(), CDkOmega, c)
      : this->F1(CDkOmega)
    );
    const volScalarField F23
    (
        kfuse ? fusedF2(k_, omega_, y_, tnu(), c) : this->F23()
    );

    // 2. Upwind weights and the convection coefficients shared by omega, k
    tmp<surfaceScalarField> tweights;
    scalarField convLower, convUpper, convDiag;

    {
        const volScalarField::Internal gamma
        (
            kfuse
          ? fusedBlend
            (
                IOobject::scopedName(this->type(), "gamma"),
                F1, gamma1_, gamma2_
            )
          : this->gamma(F1)
        );
        const volScalarField::Internal beta
        (
            kfuse
          ? fusedBlend
            (
                IOobject::scopedName(this->type(), "beta"),
                F1, beta1_, beta2_
            )
          : this->beta(F1)
        );

        if (kfuse)
        {
            fusedGbyNu(GbyNu0, omega_, F23, S2, c);
        }
        else
        {
            GbyNu0 = GbyNu(GbyNu0, F23(), S2());
        }

        tmp<volScalarField> tDomegaEff
        (
            kfuse
          ? fusedDEff
            (
                "DomegaEff", F1, alphaOmega1_, alphaOmega2_, this->nut_, tnu()
            )
          : DomegaEff(F1)
        );

        // = fvm::div(alphaRhoPhi, omega_) with bounded Gauss upwind
        tweights = pos0(alphaRhoPhi);
        tmp<fvScalarMatrix> convOmega
        (
            upwindDiv(alphaRhoPhi, tweights(), omega_)
          - fvm::Sp(tdivPhi(), omega_)
        );

        convLower = convOmega().lower();
        convUpper = convOmega().upper();
        convDiag = convOmega().diag();

        // Turbulent frequency equation
        tmp<fvScalarMatrix> omegaEqn
        (
            fvm::ddt(alpha, rho, omega_)
          + convOmega
          - fvm::laplacian(alpha*rho*tDomegaEff, omega_)
         ==
            alpha()*rho()*gamma*GbyNu0
          - fvm::SuSp((2.0/3.0)*alpha()*rho()*gamma*divU, omega_)
          - fvm::Sp(alpha()*rho()*beta*omega_(), omega_)
          - fvm::SuSp
            (
                alpha()*rho()*(F1() - scalar(1))*CDkOmega()/omega_(),
                omega_
            )
          + alpha()*rho()*beta*sqr(omegaInf_)
          + Qsas(S2(), gamma, beta)
          + omegaSource()
          + fvOptions(alpha, rho, omega_)
        );

        omegaEqn.ref().relax();
        fvOptions.constrain(omegaEqn.ref());
        omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef());
        solve(omegaEqn);
        fvOptions.correct(omega_);
        bound(omega_, this->omegaMin_);
    }

    {
        // = fvm::div(alphaRhoPhi, k_) = upwind(k) - Sp(divPhi, k): the
        // internal coefficients are omega's; subtracting the Sp matrix leaves
        // the boundary coefficients and the (zero) source bitwise unchanged
        tmp<fvScalarMatrix> convK
        (
            upwindDiv
            (
                alphaRhoPhi, tweights(), k_,
                &convLower, &convUpper, &convDiag
            )
        );
        tweights.clear();

        tmp<volScalarField> tDkEff
        (
            kfuse
          ? fusedDEff("DkEff", F1, alphaK1_, alphaK2_, this->nut_, tnu())
          : DkEff(F1)
        );

        // Turbulent kinetic energy equation
        tmp<fvScalarMatrix> kEqn
        (
            fvm::ddt(alpha, rho, k_)
          + convK
          - fvm::laplacian(alpha*rho*tDkEff, k_)
         ==
            alpha()*rho()*Pk(G)
          - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, k_)
          - fvm::Sp(alpha()*rho()*epsilonByk(F1, tgradU()), k_)
          + alpha()*rho()*betaStar_*omegaInf_*kInf_
          + kSource()
          + fvOptions(alpha, rho, k_)
        );

        tgradU.clear();

        kEqn.ref().relax();
        fvOptions.constrain(kEqn.ref());
        solve(kEqn);
        fvOptions.correct(k_);
        bound(k_, this->kMin_);
    }

    if (kfuse)
    {
        // = RASModels::kOmegaSST::correctNut(S2): kOmegaSSTBase::correctNut
        // (S2), then BasicTurbulenceModel::correctNut(), which is empty for
        // incompressible models
        {
            tmp<volScalarField> tF23(fusedF2(k_, omega_, y_, tnu(), c));
            this->nut_ =
                fusedNut(k_, omega_, tF23(), S2, this->nut_.dimensions(), c);
        }
        this->nut_.correctBoundaryConditions();
        fv::options::New(this->mesh_).correct(this->nut_);
    }
    else
    {
        correctNut(S2);
    }

    stepReuse::report
    (
        "komega",
        std::string("kOmegaSST: shared divU and k/omega convection")
      + (kfuse ? ", fused per-cell maths" : "")
    );
}

} // End namespace Foam


// ************************************************************************* //
