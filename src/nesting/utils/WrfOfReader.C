#include "nesting_utils.H"
#include "fvCFD.H"
#include "GeometricField.H"

using namespace Foam;
tmp<volVectorField> load_U(fvMesh& mesh, netCDF::NcFile& dataFile, size_t it, Foam::IOobject::writeOption opt)
{
      WrfCaseInfo wrfInfo;
      readWrfCaseInfo(&wrfInfo, dataFile);

      Info << "Retrieving U value at timestep index " << it << endl;
      // Get a variable to look how it behaves
      tmp<volVectorField> pVar;

      
      pVar = 
      (
        new volVectorField 
        (
          IOobject
          (
            "U",
            mesh.time().timeName(),
            mesh.time(),
            IOobject::NO_READ,
            opt
          ),
          mesh,
          dimVelocity,
          "zeroGradient"
          )
      );

      volVectorField& Var(pVar.ref());

      size_t Nx(wrfInfo.Ncellx), Ny(wrfInfo.Ncelly), Nz(wrfInfo.Ncellz);
      size_t N = Nz*Ny*Nx;
      float* tmp_u = new  float[N]; 
      float* tmp_un = new float[N]; 
      float* tmp_v = new  float[N]; 
      float* tmp_vn = new float[N]; 
      float* tmp_w = new  float[N]; 
      float* tmp_wn = new float[N]; 
          
      netCDF::NcVar U_wrf = dataFile.getVar("U");
      netCDF::NcVar V_wrf = dataFile.getVar("V");
      netCDF::NcVar W_wrf = dataFile.getVar("W");

      int cc = 0;
      U_wrf.getVar
      (
        std::vector<size_t>{it, 0, 0, 0},
        std::vector<size_t>{1, Nz, Ny, Nx}, 
        tmp_u
      );
      U_wrf.getVar
      (
        std::vector<size_t>{it, 0, 0, 1},
        std::vector<size_t>{1, Nz, Ny, Nx}, 
        tmp_un
      );
      V_wrf.getVar
      (
        std::vector<size_t>{it, 0, 0, 0},
        std::vector<size_t>{1, Nz, Ny, Nx}, 
        tmp_v
      );
      V_wrf.getVar
      (
        std::vector<size_t>{it, 0, 1, 0},
        std::vector<size_t>{1, Nz, Ny, Nx}, 
        tmp_vn
      );
      W_wrf.getVar
      (
        std::vector<size_t>{it, 0, 0, 0},
        std::vector<size_t>{1, Nz, Ny, Nx}, 
        tmp_w
      );
      W_wrf.getVar
      (
        std::vector<size_t>{it, 1, 0, 0},
        std::vector<size_t>{1, Nz, Ny, Nx}, 
        tmp_wn
      );
      for(size_t ibt = 0; ibt < wrfInfo.Ncellz; ibt++)
      {
        for(size_t isn = 0; isn < wrfInfo.Ncelly; isn++)
        {
          for(size_t iwe = 0; iwe < wrfInfo.Ncellx; iwe++)
          {
            cc = iwe + wrfInfo.Ncellx*isn + wrfInfo.Ncellx*wrfInfo.Ncelly*ibt;

            Var.primitiveFieldRef()[cc][0] = 0.5*(tmp_u[cc]+tmp_un[cc]);
            Var.primitiveFieldRef()[cc][1] = 0.5*(tmp_v[cc]+tmp_vn[cc]);
            Var.primitiveFieldRef()[cc][2] = 0.5*(tmp_w[cc]+tmp_wn[cc]);
          }
        }
      }
    Info << "Load U complete" << endl;
    return pVar;
}

Foam::tmp<Foam::volScalarField> load_var
(
  Foam::fvMesh& mesh, 
  netCDF::NcFile& dataFile, 
  const std::string& varname, 
  dimensionSet ds,
  size_t it,
  Foam::IOobject::writeOption opt
)
{
      WrfCaseInfo wrfInfo;
      readWrfCaseInfo(&wrfInfo, dataFile);

      Info << "Retrieving var value at timestep index " << it << endl;
      // Get a variable to look how it behaves
      tmp<volScalarField> pVar
      (
        new volScalarField 
        (
          IOobject
          (
            varname,
            mesh.time().timeName(),
            mesh.time(),
            IOobject::NO_READ,
            opt
          ),
          mesh,
          ds,
          "zeroGradient"
          )
      );

      volScalarField& Var(pVar.ref());

      size_t Nx(wrfInfo.Ncellx), Ny(wrfInfo.Ncelly), Nz(wrfInfo.Ncellz);
      size_t N = Nz*Ny*Nx;
      float* tmp = new float[N];
      netCDF::NcVar var_wrf = dataFile.getVar(varname);
      var_wrf.getVar
      (
        std::vector<size_t>{it, 0, 0, 0}, 
        std::vector<size_t>{1, Nz, Ny, Nx}, 
        tmp
      );
      int cc = 0;
      // int pcu = 0;
      // int pcv = 0;
      // int pcw = 0;
      for(size_t ibt = 0; ibt < wrfInfo.Ncellz; ibt++)
      {
        for(size_t isn = 0; isn < wrfInfo.Ncelly; isn++)
        {
          for(size_t iwe = 0; iwe < wrfInfo.Ncellx; iwe++)
          {
            cc = iwe + wrfInfo.Ncellx*isn + wrfInfo.Ncellx*wrfInfo.Ncelly*ibt;

            Var.primitiveFieldRef()[cc] = tmp[cc];
          }
        }
      }
    Info << "load var complete" << endl;
    return pVar;
}

Foam::tmp<Foam::volScalarField> load_2dvar
(
  Foam::fvMesh& mesh, 
  netCDF::NcFile& dataFile, 
  const std::string& varname, 
  dimensionSet ds,
  size_t it,
  Foam::IOobject::writeOption opt
)
{
    WrfCaseInfo wrfInfo;
    readWrfCaseInfo(&wrfInfo, dataFile);

    Info << "Retrieving var " << varname << " value at timestep index " << it << endl;
    // Get a variable to look how it behaves
    tmp<volScalarField> pVar
    (
      new volScalarField 
      (
        IOobject
        (
          varname,
          mesh.time().timeName(),
          mesh.time(),
          IOobject::NO_READ,
          opt
        ),
        mesh,
        dimensionedScalar(ds, 0)
        )
    );

    volScalarField& Var(pVar.ref());

    size_t Nx(wrfInfo.Ncellx), Ny(wrfInfo.Ncelly), Nz(wrfInfo.Ncellz);
    size_t N = Nx*Ny;
    float* tmp = new float[N];
    netCDF::NcVar var_wrf = dataFile.getVar(varname);
    std::cout << var_wrf.getDimCount() << std::endl;
    std::cout << var_wrf.getDims()[1].getSize() 
              << ", " 
              << var_wrf.getDims()[2].getSize()
              << std::endl;
    std::cout << N  << std::endl;
    var_wrf.getVar
    (
      std::vector<size_t>{it, 0, 0}, 
      std::vector<size_t>{1, Ny, Nx}, 
      tmp
    );
    label bottom_id = 0;
    for(label i=0; i < Var.boundaryField().size(); i++)
    {
      if(Var.boundaryField()[i].patch().name() == "bottom")
      {
        bottom_id = i;
        break;
      }
    }
    fvPatchScalarField& varPatch = Var.boundaryFieldRef()[bottom_id];
    
    int cc = 0, cci=0;
    // int pcu = 0;
    // int pcv = 0;
    // int pcw = 0;
      for(size_t isn = 0; isn < wrfInfo.Ncelly; isn++)
      {
        for(size_t iwe = 0; iwe < wrfInfo.Ncellx; iwe++)
        {
          cc = iwe + wrfInfo.Ncellx*isn;
          cci = isn + wrfInfo.Ncelly*iwe;

          varPatch[cci] = tmp[cc];
        }
      }
  Info << "load 2dvar complete" << endl;
  return pVar;
}
