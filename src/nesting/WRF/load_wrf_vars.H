runTime++;

Info << "Initializing variables with WRF data" << endl;
Info << "On U...." << endl;
phase.URef().primitiveFieldRef() = wrf.U();
phase.URef().correctBoundaryConditions();

Info << "On T...." << endl;
phase.thermoRef().T().primitiveFieldRef() = wrf.var("T.air");
phase.thermoRef().T().correctBoundaryConditions();

volScalarField& rho = phase.thermoRef().rho();
rho.primitiveFieldRef() = wrf.var("p") / (wrf.var("T.air") * dimensionedScalar(dimEnergy/(dimMass*dimTemperature), 287.05));

Info << "average(rho) = " << average(rho) << endl;
Info << "average(thermo.rho) = " << average(phase.thermoRef().rho()) << endl;
// p_rgh.primitiveFieldRef() = wrf.var("p_rgh");
// // p_rgh.correctBoundaryConditions();
// p.primitiveFieldRef() = wrf.var("p");
// Info << "On P...." << endl;
// phase.thermo().p().primitiveFieldRef() = wrf.interpolate<scalar>(mesh.cellCentres(), p_wrf.ref());
// fluid.correctThermo();
