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

Description
    Base class for the handling the solid trunks and branches of treeModel, providing the momentum & turbulence modification of the system. 
    
SourceFiles
    trunkPhysicsModel.C
\*---------------------------------------------------------------------------*/

#ifndef trunkPhysicsModel_H
#define trunkPhysicsModel_H

#include "fvCFD.H"
namespace Foam
{

class trunkPhysicsModel
{
    // Private Member Data
        
        //- Trunk Surface Model

        //- Internal source terms
        vectorField fU_;
        scalarField fk_;
        scalarField feps_;
        scalarField fomega_;
        tensorField fR_;

public:
    trunkPhysicsModel();
        
    virtual ~trunkPhysicsModel(){};
    
    //- Member Functions
    vectorField& fU(){return fU_;};
    scalarField& fk(){return fk_;};
    scalarField& feps(){return feps_;};
    tensorField& fR(){return fR_;};

};

}

#endif
