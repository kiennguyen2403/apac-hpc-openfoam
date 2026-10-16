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

    tmp<volTensorField> tgradU = fvc::grad(U);
    const volScalarField S2(this->S2(tgradU()));
    volScalarField::Internal GbyNu0(this->GbyNu0(tgradU(), S2));
    volScalarField::Internal G(this->GName(), nut*GbyNu0);

    // Update omega and G at the wall (see kOmegaSSTBase.C)
    omega_.boundaryFieldRef().updateCoeffs();
    omega_.boundaryFieldRef().template evaluateCoupled<coupledFvPatch>();

    const volScalarField CDkOmega
    (
        (2*alphaOmega2_)*(fvc::grad(k_) & fvc::grad(omega_))/omega_
    );

    const volScalarField F1(this->F1(CDkOmega));
    const volScalarField F23(this->F23());

    // 2. Upwind weights and the convection coefficients shared by omega, k
    tmp<surfaceScalarField> tweights;
    scalarField convLower, convUpper, convDiag;

    {
        const volScalarField::Internal gamma(this->gamma(F1));
        const volScalarField::Internal beta(this->beta(F1));

        GbyNu0 = GbyNu(GbyNu0, F23(), S2());

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
          - fvm::laplacian(alpha*rho*DomegaEff(F1), omega_)
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

        // Turbulent kinetic energy equation
        tmp<fvScalarMatrix> kEqn
        (
            fvm::ddt(alpha, rho, k_)
          + convK
          - fvm::laplacian(alpha*rho*DkEff(F1), k_)
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

    correctNut(S2);

    stepReuse::report("komega", "kOmegaSST: shared divU and k/omega convection");
}

} // End namespace Foam


// ************************************************************************* //
