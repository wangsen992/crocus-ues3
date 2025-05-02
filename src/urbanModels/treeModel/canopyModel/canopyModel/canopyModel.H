/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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
    Foam::canopyModel

Description
    Base class for the handling the physics & evapotranspiration of leaves & canopies, 
    which provides source terms to the treeModel that ultimately interacts with 
    the solver. Only to be called by the treeModel. 

    The key is on how leaves spatial domain are defined. With OpenFOAM utilities,
    cellSet can be used (or should ultimately be used), but the given the 
    potential large number of instances, there could be many surface files, or 
    many dictionaries. So this selection should provide multiple methods and 
    potentially easier interface to set up the vegetation coverage. 
    
SourceFiles
    canopyModel.C
\*---------------------------------------------------------------------------*/

#ifndef canopyModel_H
#define canopyModel_H

#include "fvCFD.H"
#include "runTimeSelectionTables.H"

#include "cellSetHashTables.H"

namespace Foam
{

// Forward declaration of treeModel
class treeModel;

class canopyModel
{
    
    // Private Member Data
    
protected:
      //- Reference to the treeModel 
      const treeModel& tree_;

      //- Dictionary 
      dictionary canopyModelDict_;



public:
    
    //- Runtime type information
    ClassName("canopyModel");
     
    // Declare run-time construtor selection table
    declareRunTimeSelectionTable
    (
      autoPtr, 
      canopyModel, 
      treeModel, 
      (
        const treeModel& tree
      ),
      (tree)
    );
    
    // Constructors
        //- Construct with components
        canopyModel
        (
          const treeModel& tree
        );

    // Selectors

        //- Return areference to the selected canopy physics model
        static autoPtr<canopyModel> New
        (
            const treeModel& tree
        );
            

    // Destructor
        virtual ~canopyModel();

    // Member Functions

        //- Return reference to the mesh
        const treeModel& tree() const;

        //- Return dictionary 
        const dictionary& dict() const;

        //- Read canopy properties dictionary
        virtual bool read();

        //- Return a list of field names that will be impacted
        virtual wordList addSupFields() const = 0;
        
        // Access to surface / geoemtry related properties

            //- cellLabels
            virtual const dictionary& surfaceModelDict() const = 0 ;

            virtual const labelHashSet& canopyCells() const = 0;
            virtual labelHashSet& canopyCells() = 0;

            virtual const dimensionedVectorCellSet& lad() const = 0;
            virtual dimensionedVectorCellSet& lad() = 0;
            virtual const dimensionedScalarCellSet& la() const = 0;
            virtual dimensionedScalarCellSet& la() = 0;
            virtual const dimensionedScalarCellSet& laLit() const = 0;
            virtual dimensionedScalarCellSet& laLit() = 0;
            virtual const dimensionedScalarCellSet& ldia() const = 0;
            virtual dimensionedScalarCellSet& ldia() = 0;

            virtual dimensionedScalarCellSet correctLaCov(vector direction)=0;
        // Momentum

            //- Return the source term to momentumTransport
            virtual dimensionedVectorCellSet& Fu() = 0;
            virtual const dimensionedVectorCellSet& Fu() const = 0;
            virtual dimensionedScalarCellSet& Fturb(const word&) = 0;
            virtual const dimensionedScalarCellSet& Fturb(const word&) const = 0;
            virtual void correctMomentumTransfer() = 0;
            //
         // Energy
            //- Return the source term to Energy
            virtual dimensionedScalarCellSet& Fhe() = 0;
            virtual const dimensionedScalarCellSet& Fhe() const = 0;
            virtual dimensionedScalarCellSet& Fq() = 0;
            virtual const dimensionedScalarCellSet& Fq() const = 0;

            //- Thermo related
            virtual dimensionedScalarCellSet& Tleaf() = 0;
            virtual const dimensionedScalarCellSet& Tleaf() const = 0;

            //- Radiation related 
            virtual dimensionedScalarCellSet& a() = 0;
            virtual const dimensionedScalarCellSet& a() const = 0;
            virtual dimensionedScalarCellSet& e() = 0;
            virtual const dimensionedScalarCellSet& e() const = 0;
            virtual dimensionedScalarCellSet& E() = 0;
            virtual const dimensionedScalarCellSet& E() const = 0;

            virtual void correctEnergyTransfer() = 0;

 
 
//         // Species
//             //- Return the source term to Yi
//             virtual tmp<volScalarField> SupYi(const word& name) = 0;
//             virtual tmp<volScalarField> SupRhoYi(const word& name) = 0;
};

}
    
#endif
