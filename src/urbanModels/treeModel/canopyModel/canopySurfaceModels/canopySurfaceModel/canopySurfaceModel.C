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
    Foam::canopyGeometryModel

\*---------------------------------------------------------------------------*/

#include "canopySurfaceModel.H"
#include "dimensionSets.H"
#include "zero.H"

namespace Foam
{

template<class BaseCanopyModel>
canopySurfaceModel<BaseCanopyModel>::canopySurfaceModel
(
    const treeModel& tree
)
:
    BaseCanopyModel(tree),
    canopyCellsIndex_(),
    lad_
    (
      IOobject
      (
        "canopy_lad",
        tree.mesh().time().constant(),
        tree.mesh(),
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
      ),
      tree.mesh(),
      dimensionedScalar(dimArea/dimVol, 0)
     ),
    la_
    (
      IOobject
      (
        "canopy_la",
        tree.mesh().time().constant(),
        tree.mesh(),
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
      ),
      tree.mesh(),
      dimensionedScalar(dimArea/dimVol, 0)
    ),
    laLit_
    (
      IOobject
      (
        "canopy_laLit",
        tree.mesh().time().constant(),
        tree.mesh(),
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
      ),
      tree.mesh(),
      dimensionedScalar(dimArea/dimVol, 0)
    ),
    laCov_
    (
      IOobject
      (
        "canopy_laCov",
        tree.mesh().time().constant(),
        tree.mesh(),
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
      ),
      tree.mesh(),
      dimensionedScalar(dimless, 0)
    ),
    ldia_
    (
      IOobject
      (
        "canopy_ldia",
        tree.mesh().time().constant(),
        tree.mesh(),
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
      ),
      tree.mesh(),
      dimensionedScalar(dimLength, 0)
    )
{
};

}
