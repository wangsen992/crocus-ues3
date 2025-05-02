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
\*---------------------------------------------------------------------------*/

#include "treeModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //
namespace Foam
{

    defineTypeNameAndDebug(treeModel, 0);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::treeModel::treeModel
(
    const fvMesh& mesh,
    const dictionary& treeModelDict
)
:
    mesh_(mesh),
    treeModelDict_(treeModelDict),
    canopy_
    (
      canopyModel::New(*this)
    ),
    trunk_
    (
      trunkModel::New(*this)
    )
    thermo_(mesh_.lookupObjectRef<fluidAtmThermo>(fluidAtmThermo::dictName))
{
}
    

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
// ************************************************************************* //
}
