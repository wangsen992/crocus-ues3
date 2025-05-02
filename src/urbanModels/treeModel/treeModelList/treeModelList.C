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

#include "IOobject.H"
#include "treeModelList.H"

namespace Foam
{

void treeModelList::init()
{
    Info << "treeModelList init…" << endl;
    if (type_ == "singleTriSurface")
    {
        this->append
        (
          new treeModel
          (
            mesh_,
            treeModelListDict_.subDict(type_ + "Coeffs")
          )
        );
    }
    Info << "treeModelList init complete." << endl;
}

treeModelList::treeModelList
(
    const fvMesh& mesh,
    const dictionary& treeModelListDict
)
:
    PtrList<treeModel>(),
    mesh_(mesh),
    treeModelListDict_(treeModelListDict),
    name_(treeModelListDict_.lookup<word>("name")),
    type_
    (
      treeModelListDict_.lookupOrDefault
      (
        "type",
        word("singleTriSurface")
      )
    ),
    lad_
    (
      IOobject
      (
        "lad",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
      ),
      mesh,
      dimensionedVector(dimArea/dimVolume, vector(0, 0, 0))
    ),
    Tleaf_
    (
      IOobject
      (
        "Tleaf",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
      ),
      mesh,
      dimensionedScalar(dimTemperature, 0)
    ),
    aTree_
    (
      IOobject
      (
        "aTree",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
      ),
      mesh,
      dimensionedScalar(dimArea/dimVolume, 0)
    ),
    eTree_
    (
      IOobject
      (
        "eTree",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
      ),
      mesh,
      dimensionedScalar(dimArea/dimVolume, 0)
    ),
    ETree_
    (
      IOobject
      (
        "ETree",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
      ),
      mesh,
      dimensionedScalar(dimMass/dimLength/pow3(dimTime), 0)
    )
{
    init();
    correct();
}

void treeModelList::correct()
{
  Info << "treeModelList correct…" << endl;
  forAllIter(PtrList<treeModel>, *this, iter)
  {
    Info << "Correcting tree 1" << endl;
    treeModel& treei(*iter);
    treei.canopy().correctMomentumTransfer();
    treei.canopy().correctEnergyTransfer();
    forAllConstIter(dimensionedVectorCellSet, treei.canopy().lad(), jter)
    {
        lad_[jter.key()] = treei.canopy().lad()[jter.key()].value();
        Tleaf_[jter.key()] = treei.canopy().Tleaf()[jter.key()].value();
        aTree_[jter.key()] = treei.canopy().a()[jter.key()].value();
        eTree_[jter.key()] = treei.canopy().e()[jter.key()].value();
        ETree_[jter.key()] = treei.canopy().E()[jter.key()].value();
    }
    Info << "Correcting tree 1 complete" << endl;
  }
  Info << "Correcting treeModelList complete." << endl;
}
}
