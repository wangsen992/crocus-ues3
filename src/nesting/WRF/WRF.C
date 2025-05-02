#include "WRF.H"
#include "treeDataFace.H"
#include "indexedOctree.H"
#include "treeBoundBox.H"

using namespace Foam;

WRF::WRF
(
  const string& ncfile_path, 
  const string& wrfCaseRoot, 
  const string& wrfCaseName, 
  const fvMesh& mesh,
  const string& foam_proj4,
  const scalar& dt,
  const scalar& T0,
  const scalar& P0
)
:
  regIOobject
  (
    IOobject
    (
      "WRF",
      wrfCaseRoot/"constant",
      mesh.time(),
      IOobject::NO_READ,
      IOobject::NO_WRITE
    )
  ),
  nc_(ncfile_path, netCDF::NcFile::read),
  runTime_(Time::controlDictName, wrfCaseRoot, wrfCaseName),
  foamMesh_(mesh),
  foamTime_(foamMesh_.time()),
  pmesh_
  (
    fvmeshFromNc(nc_, runTime_)
  ),
  ptransformer_(nullptr),
  pitransformer_(nullptr),
  dt_(dt),
  T0_(dimTemperature, T0),
  P0_(dimPressure, P0),
  thermo_
  (
    foamMesh_.lookupObjectRef<fluidAtmThermo>
    (
      IOobject::groupName
      (
        "thermophysicalProperties",
        "air"
      )
    )
  ),
  U_
  (
    IOobject
    (
      "U.air",
      runTime_.timeName(),
      runTime_,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
    ),
    pmesh_(),
    dimVelocity,
    "zeroGradient"
  )
{
  Info << "Running WRF init script" << endl;
  // load the transformation
  WRF_PROJ_PARAMS params;
  getProjAtts(nc_, params);
  std::shared_ptr<OGRSpatialReference> pcrs_wrf = getCRS(params);
  std::shared_ptr<OGRSpatialReference> pcrs_foam
  (
    new OGRSpatialReference
  );
  pcrs_foam->SetFromUserInput(foam_proj4.c_str());
  ptransformer_ = OGRCreateCoordinateTransformation(pcrs_foam.get(), pcrs_wrf.get());
  pitransformer_ = ptransformer_->GetInverse();
  Info << "transform created" << endl;
  pmesh_->movePoints
  (
    itransform(pmesh_->points())
  );
  pmesh_->setInstance(runTime_.constant());
  pmesh_->write();
  Info << "WRF init script complete (constructed & transformed to Foam CRS" << endl;
  tsearcher_.set
  (
    new meshSearch(pmesh_())
  );


  // Initialize variable resources
    volScalarFieldPtrTable_.set
    (
      "T.air",
      autoPtr<volScalarField>
      (
        new volScalarField
        (
          IOobject
          (
            "T.air",
            runTime_.timeName(),
            runTime_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
          ),
          pmesh_(),
          dimTemperature,
          "zeroGradient"
        )
      )
    );
    volScalarFieldPtrTable_.set
    (
      "H2O.air",
      autoPtr<volScalarField>
      (
        new volScalarField
        (
          IOobject
          (
            "H2O.air",
            runTime_.timeName(),
            runTime_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
          ),
          pmesh_(),
          dimless,
          "zeroGradient"
        )
      )
    );
    volScalarFieldPtrTable_.set
    (
      "thermo:rho.air",
      autoPtr<volScalarField>
      (
        new volScalarField
        (
          IOobject
          (
            "thermo.rho.air",
            runTime_.timeName(),
            runTime_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
          ),
          pmesh_(),
          dimDensity,
          "zeroGradient"
        )
      )
    );
    volScalarFieldPtrTable_.set
    (
      "p",
      autoPtr<volScalarField>
      (
        new volScalarField
        (
          IOobject
          (
            "p",
            runTime_.timeName(),
            runTime_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
          ),
          pmesh_(),
          dimPressure,
          "zeroGradient"
        )
      )
    );
}

tmp<volVectorField> WRF::read_U(size_t it, IOobject::writeOption opt)
{
  tmp<volVectorField> pu = load_U(mesh(), nc_, it, opt);
  return pu;
}

tmp<volScalarField> WRF::read_var(const word& name, dimensionSet ds, size_t it, IOobject::writeOption opt)
{
  Info << "Start loading var" << endl;
  tmp<volScalarField> pvar = load_var(mesh(), nc_, name, ds, it, opt);
  pvar->correctBoundaryConditions();
  return pvar;
}

tmp<volScalarField> WRF::read_var2d(const word& name, dimensionSet ds, size_t it, IOobject::writeOption opt)
{
  Info << "Start loading var" << endl;
  tmp<volScalarField> pvar = load_2dvar(mesh(), nc_, name, ds, it, opt);
  return pvar;
}

void Foam::WRF::updateVars(label it)
{
    runTime_.setTime(this->dt() * it, it);
    U_ = this->read_U(it);
    U_.correctBoundaryConditions();
    volScalarFieldPtrTable_["p"]() = this->read_var("P", dimPressure, it) + this->read_var("PB", dimPressure, it);
    volScalarFieldPtrTable_["p"]().correctBoundaryConditions();

    volScalarFieldPtrTable_["T.air"]() 
      = (this->read_var("T", dimTemperature, it) + this->T0()) 
        * pow
          (
            volScalarFieldPtrTable_["p"]()/this->P0(), 
            0.286
          );
    volScalarFieldPtrTable_["T.air"]().correctBoundaryConditions();
    volScalarFieldPtrTable_["H2O.air"]() = this->read_var("QVAPOR", dimless, it);
    volScalarFieldPtrTable_["H2O.air"]().correctBoundaryConditions();
    volScalarFieldPtrTable_["thermo:rho.air"]() = volScalarFieldPtrTable_["p"]() / (volScalarFieldPtrTable_["T.air"]() * dimensionedScalar(dimEnergy/(dimMass*dimTemperature), 287.05));
    volScalarFieldPtrTable_["thermo:rho.air"]().correctBoundaryConditions();

    // Debug for WRF data
    // runTime_++;
    // runTime_.write();

}

point WRF::transform(const point& pt)
{
  double x(pt.x()), y(pt.y()), z(pt.z());
  ptransformer_->Transform(1, &x, &y);
  return point{x,y,z};
}

pointField WRF::transform(const pointField& pts)
{
  pointField outPts(pts.size());
  for(int i=0; i<pts.size(); i++)
  {
    outPts[i] = transform(pts[i]);
  }
  return outPts;
}

point WRF::itransform(const point& pt)
{
  double x(pt.x()), y(pt.y()), z(pt.z());
  pitransformer_->Transform(1, &x, &y);
  return point{x,y,z};
}


pointField WRF::itransform(const pointField& pts)
{
  pointField outPts(pts.size());
  for(int i=0; i<pts.size(); i++)
  {
    outPts[i] = itransform(pts[i]);
  }
  return outPts;
}

void WRF::terraform_to_wrf(fvMesh& mesh)
{
  Info << "Terraforming to WRF" << endl;
  fvMesh& wrfMesh(this->mesh());
  const polyPatch& wrfGroundPatch(wrfMesh.boundaryMesh()["bottom"]);

  indexedOctree<treeDataFace> wrfTree
  (
    treeDataFace(false, wrfGroundPatch),
    treeBoundBox(boundBox(wrfGroundPatch.localPoints())),
    10,
    10,
    3
  );

  const polyPatch& foamGroundPatch(mesh.boundaryMesh()["bottom"]);
  indexedOctree<treeDataFace> foamTree
  (
    treeDataFace(false, foamGroundPatch),
    treeBoundBox(boundBox(foamGroundPatch.points())),
    10,
    10,
    3
  );
  // Interp the ground points from openfoam to wrf bottom patch
  Foam::vector ll{0,0,10000};
  Foam::vector zvec{0,0,1};
  auto findVec = [&](const point& pt)
  {
    auto wrfHit =  wrfTree.findLine(pt-ll, pt+ll);
    if (!wrfHit.hit())
    {
      Info << "WRF bounds " << wrfTree.bb() << endl;
      Info << "WRF not hit " << pt << endl;
    }
    point wrfPt = wrfHit.hitPoint(); 
    auto foamHit = foamTree.findLine(pt-ll, pt+ll);
    point foamPt;
    if (!foamHit.hit())
    {
      foamPt = foamTree.findNearest(pt, 1e9).hitPoint();
    }
    else
    {
      foamPt = foamHit.hitPoint(); 
    }

    return zvec * wrfPt.z();
    // compress the high altitude regions
  };

  {
    pointField foamPts = mesh.points();
    vectorField vec(foamPts.size());
    std::transform
    (
      foamPts.cbegin(), 
      foamPts.cend(), 
      vec.begin(),
      findVec
    );

    // assuming foam mesh is always lower than wrf mesh
    scalar zmax = gMax(foamPts.component(2));
    scalar zmin = gMin(foamPts.component(2));
    scalar vec_zmax(gMax(vec.component(2)));
    scalar vec_zmin(gMin(vec.component(2)));
    Info << vec_zmax - vec_zmin << endl;

    vec = Foam::vector{0,0,vec_zmin} 
        +(
            (zmax - foamPts.component(2))/(zmax-zmin)
           *(vec - Foam::vector{0,0,vec_zmin})
         );

    mesh.movePoints(foamPts + vec); // Some points are
                                                      // outside the wrf domain
    mesh.setInstance(mesh.time().constant());
    mesh.write();
    mesh.moving(false); // set to false so solver doesn't require V0
  Info << "Terraforming to WRF complete" << endl;
  }
}
