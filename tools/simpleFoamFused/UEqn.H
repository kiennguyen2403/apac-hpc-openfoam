    // Momentum predictor

    MRF.correctBoundaryVelocity(U);

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::div(phi, U)
      + MRF.DDt(U)
      + turbulence->divDevReff(U)
     ==
        fvOptions(U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (simple.momentumPredictor())
    {
        if (!simpleFoamFused::off())
        {
            // = solve(UEqn == -fvc::grad(p)) without copying the matrix: the
            // source with V*(-grad(p)) is swapped in for the solve
            simpleFoamFused::report("momentum solve without matrix copy");

            vectorField b
            (
                simpleFoamFused::momentumSource(UEqn, fvc::grad(p)())
            );
            UEqn.source().swap(b);
            UEqn.solve();
            UEqn.source().swap(b);
        }
        else
        {
            solve(UEqn == -fvc::grad(p));
        }

        fvOptions.correct(U);
    }
