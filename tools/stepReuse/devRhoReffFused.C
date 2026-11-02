/*---------------------------------------------------------------------------*\
    Replacement for linearViscousStress<incompressible RAS>::devRhoReff(U)
    (explicit specialisation; reached through the vtable under LD_PRELOAD).
    Called once per time step by the wallShearStress function object
    (devReff() -> devRhoReff() -> devRhoReff(U)); forces/forceCoeffs use
    their own patch-only evaluation.

    Stock:  (-(alpha*rho*nuEff()))*devTwoSymm(fvc::grad(U))
    with nuEff() = nut() + nu() (RASModel) and alpha, rho one-fields: four
    full-mesh field operations (nuEff, devTwoSymm, negation, product), each
    allocating a cell field and writing it back (about 420-520 B/cell, job
    182232191: 1.4% of the step at 4 nodes).

    Here one cell loop and one loop per patch evaluate
        (-(nut[i] + nu[i]))*devTwoSymm(gradU[i])
    with the same inline per-element functions in the same order, into the
    result field directly (about 180 B/cell): bit-identical. grad(U) comes
    from fvc::grad(U) as in the stock code (the shipped fvSolution caches it).
    The function object itself is untouched: it runs every step and gets the
    same field.

    STEP_REUSE_OFF=devrho runs the stock function.
\*---------------------------------------------------------------------------*/

#include "IncompressibleTurbulenceModel.H"
#include "incompressible/transportModel/transportModel.H"
#include "RASModel.H"
#include "linearViscousStress.H"
#include "fvc.H"
#include "stepReuse.H"

namespace Foam
{

typedef RASModel<IncompressibleTurbulenceModel<transportModel>>
    stepReuseIncompressibleRAS;

template<>
tmp<volSymmTensorField>
linearViscousStress<stepReuseIncompressibleRAS>::devRhoReff
(
    const volVectorField& U
) const
{
    if (stepReuse::off("devrho"))
    {
        stepReuse::report("devrho", "devRhoReff: stock");

        typedef tmp<volSymmTensorField> (*devRhoReffFn)
        (
            const linearViscousStress<stepReuseIncompressibleRAS>*,
            const volVectorField&
        );
        static const devRhoReffFn stock = reinterpret_cast<devRhoReffFn>
        (
            stepReuse::next
            (
                "_ZNK4Foam19linearViscousStressINS_8RASModelINS_"
                "29IncompressibleTurbulenceModelINS_14transportModelEEEEEE"
                "10devRhoReffERKNS_14GeometricFieldINS_6VectorIdEENS_"
                "12fvPatchFieldENS_7volMeshEEE"
            )
        );
        return stock(this, U);
    }

    stepReuse::report("devrho", "devRhoReff: one-pass evaluation");

    tmp<volScalarField> tnut(this->nut());
    tmp<volScalarField> tnu(this->nu());
    const volScalarField& nut = tnut();
    const volScalarField& nu = tnu();

    tmp<volTensorField> tgradU(fvc::grad(U));
    const volTensorField& gradU = tgradU();

    tmp<volSymmTensorField> tres
    (
        new volSymmTensorField
        (
            IOobject
            (
                IOobject::groupName("devRhoReff", this->alphaRhoPhi_.group()),
                this->runTime_.timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh_,
            dimensioned<symmTensor>
            (
                nut.dimensions()*gradU.dimensions(),
                Zero
            )
        )
    );
    volSymmTensorField& res = tres.ref();

    // Cells: nuEff = nut + nu, then (-nuEff)*devTwoSymm(gradU), as the
    // stock field operations do per element
    {
        const scalarField& nuti = nut.primitiveField();
        const scalarField& nui = nu.primitiveField();
        const tensorField& G = gradU.primitiveField();
        symmTensorField& r = res.primitiveFieldRef();

        forAll(r, celli)
        {
            r[celli] = (-(nuti[celli] + nui[celli]))*devTwoSymm(G[celli]);
        }
    }

    // Patches (calculated), from the patch values of nut, nu and grad(U)
    volSymmTensorField::Boundary& rb = res.boundaryFieldRef();

    forAll(rb, patchi)
    {
        const scalarField& nutp = nut.boundaryField()[patchi];
        const scalarField& nup = nu.boundaryField()[patchi];
        const tensorField& Gp = gradU.boundaryField()[patchi];
        symmTensorField& rp = rb[patchi];

        forAll(rp, i)
        {
            rp[i] = (-(nutp[i] + nup[i]))*devTwoSymm(Gp[i]);
        }
    }

    return tres;
}

} // End namespace Foam


// ************************************************************************* //
