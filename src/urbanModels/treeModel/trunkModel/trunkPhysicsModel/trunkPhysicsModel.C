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

#include "trunkPhysicsModel.H"

namespace Foam
{

trunkPhysicsModel::trunkPhysicsModel()
:
fU_(0),
fk_(0),
feps_(0),
fomega_(0),
fR_(0)
{}

}
