{
    alphat = turbulence->nut()/Prt;
    alphat.correctBoundaryConditions();

    volScalarField alphaEff("alphaEff", turbulence->nu()/Pr + alphat);

    fvScalarMatrix qEqn
    (
        fvm::ddt(q)
      + fvm::div(phi, q)
      - fvm::laplacian(alphaEff, q)
     ==
        fvOptions(q)
      + Source_dqtdt_advec
      + tau * (Source_qt_LS - q)
    );

    qEqn.relax();

    fvOptions.constrain(qEqn);

    qEqn.solve();

    fvOptions.correct(q);
}
