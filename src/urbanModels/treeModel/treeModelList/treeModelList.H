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
    Foam::treeModelList

Description
    Class for organizing all trees into a bulk form
    
SourceFiles
    treeModelList.C
\*---------------------------------------------------------------------------*/
#ifndef treeModelList_H
#define treeModelList_H

#include "volFieldsFwd.H"
#include "PtrList.H"
#include "treeModel.H"

namespace Foam
{

class treeModelList
  : public PtrList<treeModel>
{
    // Private Data

        //- Reference to mesh
        const fvMesh& mesh_;

        // Dictionary that holds all information of this treeModelList
        dictionary treeModelListDict_;

        // name of the treeModelList
        word name_;

        // type
        word type_;

        // Tree-relevant volField variables that should be registered.

            // Leaf are density vector
            volVectorField  lad_;

            // Leaf are density vector
            volScalarField  Tleaf_;

            // absorption coefficient from trees
            volScalarField  aTree_;
        
            // emission coefficient from trees
            volScalarField  eTree_;
        
            // emission contribution from trees
            volScalarField  ETree_;
    // Private Member Functions

          void init();
  
public:

    // Constructor

        // Construct from mesh directly with dictionary and mesh
        treeModelList
        (
          const fvMesh& mesh,
          const dictionary& treeModelListDict
        );

    // Destructor
    virtual ~treeModelList(){}

    // Member Function

      //- Access

        word name(){return name_;};
        word treeType(){return type_;};
        dictionary dict(){return treeModelListDict_;};

        volVectorField& lad(){return lad_;};
        const volVectorField& lad() const {return lad_;};

        volScalarField& a(){return aTree_;};
        const volScalarField& a() const {return aTree_;};

        volScalarField& e(){return eTree_;};
        const volScalarField& e() const {return eTree_;};

        volScalarField& E(){return ETree_;};
        const volScalarField& E() const {return ETree_;};

        volScalarField& Tleaf(){return Tleaf_;};
        const volScalarField& Tleaf() const {return Tleaf_;};
        
      //- 
        void correct();
};
}

#endif
