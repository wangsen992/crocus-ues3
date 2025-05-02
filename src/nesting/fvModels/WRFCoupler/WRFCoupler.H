/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2021 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::fv::WRFCoupler

Description
    WRF coupler to be used for the atm*Foam solvers loading input 

    Weight for each cell is a summation of nudging and lateral boundary coefficients

Usage
    Example usage:
    \verbatim
    buoyancyForce1
    {
        type        WRFCoupler;

        U           U;            // Name of the velocity field
        Ug          (10 0 0);     // Value of geostrophic velocity
        f           1e-4;         // Value of Coriolis Parameter
    }
    \endverbatim

SourceFiles
    WRFCoupler.C

\*---------------------------------------------------------------------------*/

#ifndef WRFCoupler_H
#define WRFCoupler_H

#include "HashTable.H"
#include "fvModel.H"
#include "uniformDimensionedFields.H"
#include "dimensionedVector.H"
#include "fluidAtmThermo.H"
#include "nesting_utils.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                        Class buoyancyForce Declaration
\*---------------------------------------------------------------------------*/
template<typename Type>
tmp<GeometricField<Type, fvPatchField, volMesh>> smooth(GeometricField<Type, fvPatchField, volMesh>& psi, label smooth_n=5);
 
class WRFCoupler
:
    public fvModel
{
    // Private Data

        //- Host case fvMesh
        const fvMesh& mesh_;
        Time wrfTime_;
        autoPtr<fvMesh> pwrfMesh_;
        const fluidAtmThermo& thermo_;

        dictionary varList_;

        //- CellList for applying nesting
        DynamicList<label, 20> nestingCells_;
        HashTable<scalar, label> nestingCellTbl_;
        pointField nestingCellCentres_;
        volScalarField cellWeights_;
        scalar nestingDist_;
        scalar nestingDistTop_;
        scalar relaxationFactor_;

        //- Flags to control update
        label currTimeInd_;

        //- Optional phase name
        word phaseName_;

        //- 
        label wrfi_;
        volVectorField projUold_;
        volVectorField projUnew_;
        HashTable<autoPtr<volScalarField>, word> projVolScalarFieldPtrTableOld_;
        HashTable<autoPtr<volScalarField>, word> projVolScalarFieldPtrTableNew_;

    // Private Member Functions

        //- Non-virtual read
        void readCoeffs();

        //- static function for smoothing

public:

    //- Runtime type information
    TypeName("WRFCoupler");

    // Constructors

        //- Construct from explicit source name and mesh
        WRFCoupler
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Disallow default bitwise copy construction
        WRFCoupler(const WRFCoupler&) = delete;

    // Member Functions

        // Checks

            //- Return the list of fields for which the fvModel adds source term
            //  to the transport equation
            virtual wordList addSupFields() const;


        // Evaluate


            //- Add explicit contribution to phase momentum equation
            template<typename T>
            void addSup
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                fvMatrix<T>& eqn,
                const word& fieldName
            ) const;
            
           // interpolate multiple points from source field
           template<typename Type>
           static Field<Type> interpolate(const Field<point>& pts, const GeometricField<Type, fvPatchField, volMesh>& psi, const word& interpMethod="cell");
 
        
        // Correct fvModel
        void correct();


        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const WRFCoupler&) = delete;
        
        void addSup
        (
            const volScalarField& rho,
            fvMatrix<Foam::vector>& eqn,
            const word& fieldName
        ) const;
        void addSup
        (
            const volScalarField& rho,
            fvMatrix<Foam::scalar>& eqn,
            const word& fieldName
        ) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<typename Type>
Field<Type> fv::WRFCoupler::interpolate(const Field<point>& pts, const GeometricField<Type, fvPatchField, volMesh>& psi, const word& interpMethod)
{
  meshSearch searcher(psi.mesh());

  autoPtr<interpolation<Type>> interp
  (
    interpolation<Type>::New(interpMethod, psi)
  );

  pointField wrf_pts = pts;
  Field<Type> interpVals(pts.size(), Zero);
  Info << "[WRF] gAverage(interpVals) " << gAverage(interpVals)
       << "; gAverage(interp->psi) = " << gAverage(interp->psi()) << endl;
  for(size_t i=0; i < pts.size(); i++)
  {
    label cellInd = searcher.findCell(wrf_pts[i]);
    if(cellInd < 0)
    {
      Info << "[Debug] point not in domain: " << wrf_pts[i];
      cellInd = searcher.findNearestCell(wrf_pts[i]);
      Info << " distance = " << searcher.mesh().points()[cellInd].z() - wrf_pts[i].z() << endl;
    }
    interpVals[i] = interp->interpolate
    (
      wrf_pts[i],
      cellInd
    );
  }
  interp.clear();
  Info << "[WRF] gAverage(psi) = " << gAverage(psi) 
       << "gAverage(interpVals) = " << gAverage(interpVals)
       << endl;
  return interpVals;
}

template<typename Type>
tmp<GeometricField<Type, fvPatchField, volMesh>> smooth(GeometricField<Type, fvPatchField, volMesh>& psi, label smooth_n)
{
    typedef GeometricField<Type, fvPatchField, volMesh> psiType;
    psiType f0(psi);
    psiType f(psi);
    label smoothingPasses = smooth_n;
    for (label pass=1; pass<=smoothingPasses; pass++)
    {
        f = fvc::average(fvc::interpolate(f0));
        f.correctBoundaryConditions();
        f0 = f;
    }
    return tmp<psiType> ( new psiType(f0));
}

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

