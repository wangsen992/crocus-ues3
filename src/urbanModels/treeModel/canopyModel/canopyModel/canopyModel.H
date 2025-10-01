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

#include "volFieldsFwd.H"
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

            virtual const labelHashSet& canopyCells() const = 0;
            virtual labelHashSet& canopyCells() = 0;

            virtual const volScalarField& lad() const = 0;
            virtual volScalarField& lad() = 0;
            virtual const volScalarField& la() const = 0;
            virtual volScalarField& la() = 0;
            virtual const volScalarField& laLit() const = 0;
            virtual volScalarField& laLit() = 0;
            virtual const volScalarField& ldia() const = 0;
            virtual volScalarField& ldia() = 0;

            virtual volScalarField correctLaCov(vector direction)=0;
        // Momentum

            //- Return the source term to momentumTransport
            virtual volVectorField& Fu() = 0;
            virtual const volVectorField& Fu() const = 0;
            virtual volScalarField& Fturb(const word&) = 0;
            virtual const volScalarField& Fturb(const word&) const = 0;
            virtual void correctMomentumTransfer() = 0;
            //
         // Energy
            //- Return the source term to Energy
            virtual volScalarField& FT() = 0;
            virtual const volScalarField& FT() const = 0;
            virtual volScalarField& Fq() = 0;
            virtual const volScalarField& Fq() const = 0;

            //- Thermo related
            virtual volScalarField& Tleaf() = 0;
            virtual const volScalarField& Tleaf() const = 0;

            //- Radiation related 
            virtual volScalarField& a() = 0;
            virtual const volScalarField& a() const = 0;
            virtual volScalarField& e() = 0;
            virtual const volScalarField& e() const = 0;
            virtual volScalarField& E() = 0;
            virtual const volScalarField& E() const = 0;

            virtual void correctEnergyTransfer() = 0;

 
 
//         // Species
//             //- Return the source term to Yi
//             virtual tmp<volScalarField> SupYi(const word& name) = 0;
//             virtual tmp<volScalarField> SupRhoYi(const word& name) = 0;
};

}
    
#endif
