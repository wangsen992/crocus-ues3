// Assuming foam mesh is already created
Info << "Creating WRF object..." << endl;
IOdictionary wrf_dict
(
  IOobject
  (
    "WRFDict",
    runTime.constant(),
    runTime,
    IOobject::MUST_READ
  )
);
Info << wrf_dict << endl;
WRF wrf
(
  wrf_dict["ncfile"], 
  wrf_dict["wrfRoot"], 
  "wrfCase",
  mesh,
  wrf_dict["proj4"],
  wrf_dict.lookupOrAddDefault("dt", 3600, false, false),
  wrf_dict.lookupOrAddDefault("T0", 290, false, false),
  wrf_dict.lookupOrAddDefault("P0", 1e5, false, false)
);

// When in parallel, terraform leads to mismatch at the proc boundaries
// by a small amount as the same point is translated by 
// slightly different amount
wrf.terraform_to_wrf(mesh);
wrf.updateVars(0);

// Load variables (besides U) to ptrTable
List<dictionary> variableList
(
  wrf_dict["variables"]
);
Info << "loaded variableList" << endl;
for(dictionary vd : variableList)
{
  Info << vd.lookupOrDefault<word>("name", "null") << endl;
}

// initialize variables
phaseModel& phase = fluid.movingPhases()[0];

Info << "Initializing variables with WRF data" << endl;

// Variables used from WRF: U, H2O.air, T.air, p, e.air, rho
Info << "On U...." << endl;
phase.URef().primitiveFieldRef() = wrf.interpolate(mesh.cellCentres(), wrf.U());
phase.URef().correctBoundaryConditions();
phase.phiRef() = fvc::flux(phase.URef());
phase.alphaPhiRef() = fvc::flux(phase.URef());

Info << "On T...." << endl;
phase.thermoRef().T().primitiveFieldRef() = wrf.interpolate(mesh.cellCentres(), wrf.var("T.air"));
phase.thermoRef().T().correctBoundaryConditions();

phase.YRef()[1].primitiveFieldRef() = wrf.interpolate(mesh.cellCentres(), wrf.var("H2O.air"));
phase.YRef()[0] = dimensionedScalar(dimless, 1) - phase.YRef()[1];
phase.YRef()[0].correctBoundaryConditions();
phase.YRef()[1].correctBoundaryConditions();

// volScalarField& wrf_e
// (
//   const_cast<volScalarField&>(wrf.interpolate(mesh.cellCentres(), wrf.var("e.air")))
// );
// wrf_e = phase.thermoRef().he(phase.thermoRef().p(), phase.thermoRef().T());
phase.thermoRef().he().primitiveFieldRef() = phase.thermoRef().he(phase.thermoRef().p(), phase.thermoRef().T());
phase.thermoRef().he().correctBoundaryConditions();
