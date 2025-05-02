// WRF class as a connector from OpenFOAM to read netcdf read from WRF
// WRF reads an foamcase created with the utility foamWRF
// Used as a resource manager, and all data are stored in the master node,
// whereas the rest simply request data from the master. 

#ifndef WRF_H
#define WRF_H

#include "HashTable.H"
#include "fvMesh.H"
#include "meshSearch.H"
#include "nesting_utils.H"
#include "fluidAtmThermo.H"

namespace Foam
{
  class WRF: public regIOobject
  {
    // Private member attributes
      //- data reference
      netCDF::NcFile nc_;
      Time runTime_; // Time registry for the main of case
      const fvMesh& foamMesh_;
      const Time& foamTime_;
      autoPtr<fvMesh> pmesh_;

      //- reference 
      autoPtr<meshSearch> tsearcher_;
      OGRCoordinateTransformation* ptransformer_;
      OGRCoordinateTransformation* pitransformer_;
      scalar  dt_;
      dimensionedScalar  T0_;
      dimensionedScalar  P0_;

      //- RAII for WRF data
        //- Thermo class
        const fluidAtmThermo& thermo_;

        //- List of varialbes from netcdf
        volVectorField U_;
        HashTable<autoPtr<volScalarField>, word> volScalarFieldPtrTable_;
        
        //- Flags to control update
        label currTimeInd_;

  public:

    WRF
    (
      const string& ncfile_path, 
      const string& wrfCaseRoot, 
      const string& wrfCaseName, 
      const fvMesh& mesh,
      const string& foam_proj4,
      const scalar& dt,
      const scalar& T0,
      const scalar& P0
    );

    ~WRF(){};

    // regIOobject inheritance
      
      bool writeData(Ostream&) const
      {
        return true;
      }

    // Member functions for access

      const scalar& dt(){return dt_;};
      const dimensionedScalar& T0(){return T0_;};
      const dimensionedScalar& P0(){return P0_;};
      
      netCDF::NcFile& file(){return nc_;};

      autoPtr<fvMesh> pmesh(){return pmesh_;};

      fvMesh& mesh(){return pmesh_();};

      Time& time(){return runTime_;};

      const Time& foamTime(){return foamTime_;};

      template<typename Type>
      inline bool hasVar(const string& name)
      {
        return pmesh_->foundObject<Type>(name);
      };

    // Load variables from netcdf file to mesh

      tmp<volVectorField> read_U(size_t it, IOobject::writeOption opt=IOobject::NO_WRITE);

      tmp<volScalarField> read_var(const word& name, dimensionSet ds, size_t it, IOobject::writeOption opt=IOobject::NO_WRITE);

      tmp<volScalarField> read_var2d(const word& name, dimensionSet ds, size_t it, IOobject::writeOption opt=IOobject::NO_WRITE);

      const volVectorField& U(){return U_;}

      const volScalarField& var(const word& name)
      {
        return volScalarFieldPtrTable_[name]();
      }
      
      //- Non-virtual update variables
      void updateVars(label it);

    // Coordinate transformation
      // Transform one point from input to WRF
      point transform(const point&);

      // Transform pointField from input to WRF
      pointField transform(const pointField&);
      
      // Transform one point from input to to foam from WRF
      point itransform(const point&);

      // Transform pointField from input to foam from WRF
      pointField itransform(const pointField&);

      // terraform the openfoam case to wrf terrain
      void terraform_to_wrf(fvMesh& mesh);

    // Interpolation functions
    // interpolate one point from source field
       template<typename Type>
       Type interpolate(const point& pt, const GeometricField<Type, fvPatchField, volMesh>& psi, const word& interpMethod="pointMVC");
 
       // interpolate multiple points from source field
       template<typename Type>
       Field<Type> interpolate(const Field<point>& pts, const GeometricField<Type, fvPatchField, volMesh>& psi, const word& interpMethod="pointMVC");
 
   };
}
#include "WRFI.C"

#endif
