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
    Foam::trunkModel
\*---------------------------------------------------------------------------*/

#ifndef trunkModel_H
#define trunkModel_H

#include "fvCFD.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

class treeModel;

class trunkModel
{
    // Private Member Data

      //- Reference to Tree
      const treeModel& tree_;

      //- Dictionary
      dictionary trunkModelDict_;

public:
    
    //- Runtime type information
    ClassName("trunkModel");

    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
      autoPtr, 
      trunkModel,
      treeModel,
      (
        const treeModel& tree
      ),
      (tree)
    );

    // Constructors
    trunkModel
    (
      const treeModel& tree
    );
    
    // Selectors

        //- Return areference to the selected canopy physics model
        static autoPtr<trunkModel> New
        (
            const treeModel& tree
        );

    // Destructor
    virtual ~trunkModel();
    
    //- Return reference to the mesh
    const treeModel& tree() const;

    //- Return dictionary 
    const dictionary& dict() const;


};
}
#endif
