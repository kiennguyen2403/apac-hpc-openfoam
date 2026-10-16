{
    volScalarField rAU(1.0/UEqn.A());
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));

    tmp<volScalarField> rAtU(rAU);
    autoPtr<surfaceScalarField> phiHbyAPtr;
    bool fused = false;

    if (!simpleFoamFused::off() && simple.consistent())
    {
        rAtU = 1.0/(1.0/rAU - UEqn.H1());

        const volScalarField d(rAtU() - rAU);

        fused = simpleFoamFused::pEqnOk
        (
            mesh, p, MRF.MRFZoneList::size(), simple.consistent(), HbyA, d
        );

        if (fused)
        {
            simpleFoamFused::report("pressure equation steps fused");

            // = fvc::flux(HbyA), MRF.makeRelative, adjustPhi (no-ops here),
            //   += fvc::interpolate(rAtU() - rAU)*fvc::snGrad(p)*mesh.magSf()
            phiHbyAPtr.reset
            (
                new surfaceScalarField
                (
                    "phiHbyA",
                    simpleFoamFused::phiHbyA(HbyA, d, p)
                )
            );

            // = HbyA -= (rAU - rAtU())*fvc::grad(p)
            simpleFoamFused::HbyAUpdate(HbyA, rAU, rAtU(), fvc::grad(p)());
        }
        else
        {
            // rAtU is recomputed below, in the stock order
            rAtU = rAU;
        }
    }

    if (!fused)
    {
        phiHbyAPtr.reset(new surfaceScalarField("phiHbyA", fvc::flux(HbyA)));
        surfaceScalarField& phiHbyA = phiHbyAPtr.ref();
        MRF.makeRelative(phiHbyA);
        adjustPhi(phiHbyA, U, p);

        if (simple.consistent())
        {
            rAtU = 1.0/(1.0/rAU - UEqn.H1());
            phiHbyA +=
                fvc::interpolate(rAtU() - rAU)*fvc::snGrad(p)*mesh.magSf();
            HbyA -= (rAU - rAtU())*fvc::grad(p);
        }
    }

    surfaceScalarField& phiHbyA = phiHbyAPtr.ref();

    tUEqn.clear();

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p, U, phiHbyA, rAtU(), MRF);

    // Non-orthogonal pressure corrector loop
    while (simple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rAtU(), p) == fvc::div(phiHbyA)
        );

        pEqn.setReference(pRefCell, pRefValue);

        pEqn.solve();

        if (simple.finalNonOrthogonalIter())
        {
            if (fused && pEqn.nMatrices() == 1)
            {
                // = phiHbyA - pEqn.flux()
                phi = simpleFoamFused::phiCorrected(phiHbyA, pEqn, p);
            }
            else
            {
                phi = phiHbyA - pEqn.flux();
            }
        }
    }

    #include "continuityErrs.H"

    // Explicitly relax pressure for momentum corrector
    p.relax();

    // Momentum corrector
    if (fused)
    {
        // = HbyA - rAtU()*fvc::grad(p)
        U = simpleFoamFused::correctedVelocity(HbyA, rAtU(), fvc::grad(p)());
    }
    else
    {
        U = HbyA - rAtU()*fvc::grad(p);
    }
    U.correctBoundaryConditions();
    fvOptions.correct(U);
}
