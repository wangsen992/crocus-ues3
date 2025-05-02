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
    Foam::treeModel

Description
    A base class for vegetation related models based on fvModel, such that 
    IO and interface with the equations can be handled through fvModels in 
    a centralized manner
    
SourceFiles
    treeModel.C
\*---------------------------------------------------------------------------*/
#ifndef treeModel_H
#define treeModel_H

#include "fvCFD.H"
#include "canopyModel.H"
#include "trunkModel.H"
#include "fluidAtmThermo.H"

namespace Foam
{
class trunkModel;

class treeModel
{
    
    // Private Member Data

        //- Reference to fvMesh
        const fvMesh& mesh_;

        //- Dictionary
        const dictionary treeModelDict_;
    
        //- Pointers to components of trees
        autoPtr<canopyModel> canopy_;
        autoPtr<trunkModel>  trunk_;

        //- Reference to embedded models
        fluidAtmThermo& thermo_;

    // Private Member Functions

public: 
    
    //- Runtime type information
    TypeName("treeModel");

    // Constructors
    
        //- Constructor from components
        treeModel
        (
          const fvMesh& mesh,
          const dictionary& treeModelDict
        );

        static autoPtr<treeModel> New
        (
          const fvMesh& mesh,
          const dictionary& treeModelDict
        );
        

    //- Destructor
    virtual ~treeModel()
    {}

    //- Member Functions
    wordList addSupFields();
    const fvMesh& mesh() const {return mesh_;}
    const dictionary& dict() const {return treeModelDict_;}

    canopyModel& canopy() {return canopy_();}
    const canopyModel& canopy() const {return canopy_();}

    trunkModel& trunk() {return trunk_();}
    const trunkModel& trunk() const {return trunk_();}

    fluidAtmThermo& thermo() {return thermo_;}
    const fluidAtmThermo& thermo() const {return thermo_;}
};

}
#endif 
