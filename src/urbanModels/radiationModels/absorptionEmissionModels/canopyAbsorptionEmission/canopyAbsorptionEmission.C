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

#include "canopyAbsorptionEmission.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace radiationModels
{
namespace absorptionEmissionModels
{

    defineTypeNameAndDebug(canopy, 0);

    addToRunTimeSelectionTable
    (
        absorptionEmissionModel,
        canopy,
        dictionary
    );

canopy::canopy
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    absorptionEmissionModel(dict, mesh),
    coeffsDict_(dict.subDict(typeName + "Coeffs"))
{}


canopy::~canopy(){}

tmp<volScalarField> canopy::aDisp(const label bandI) const
{
    tmp<volScalarField> ta
    (
      volScalarField::New
      (
        "a",
        mesh_,
        dimensionedScalar(dimless/dimLength, 0)
      )
    );
    
    if (bandI == 0)
    {
      ta.ref() = mesh().lookupObjectRef<volScalarField>("aTree");
    }

    return ta;
}

tmp<volScalarField> canopy::eDisp(const label bandI) const
{
    tmp<volScalarField> te
    (
      volScalarField::New
      (
        "e",
        mesh_,
        dimensionedScalar(dimless/dimLength, 0)
      )
    );

    if (bandI == 1)
    {
      te.ref() = mesh().lookupObjectRef<volScalarField>("eTree");
    }

    return te;
}

tmp<volScalarField> canopy::EDisp(const label bandI) const
{
    tmp<volScalarField> tE
    (
      volScalarField::New
      (
        "E",
        mesh_,
        dimensionedScalar(dimMass/dimLength/pow3(dimTime), 0)
      )
    );

    if (bandI == 1)
    {
      tE.ref() = mesh().lookupObjectRef<volScalarField>("ETree");
    }
    return tE;
}


}
}
}
