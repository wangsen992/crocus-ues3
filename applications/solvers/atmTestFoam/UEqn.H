    // Solve the momentum equation

    MRF.correctBoundaryVelocity(U);

    // Info << "Before UEqn" << endl;
    fvVectorMatrix UEqn
    (
        fvm::ddt(U) + fvm::div(phi, U)
      + MRF.DDt(U)
      + turbulence->divDevReff(U)
     ==
        fvOptions(U)
      + SourceU
    );

    UEqn.relax();

    // Info << "before fvoption constraint" << endl;
    fvOptions.constrain(UEqn);
    // Info << "after fvoption constraint" << endl;

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                  - ghf*fvc::snGrad(rhok)
                  - fvc::snGrad(p_rgh)
                )*mesh.magSf()
            )
        );

        fvOptions.correct(U);
    }
