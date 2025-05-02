/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "binaryWideBand.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{
namespace absorptionEmissionModels
{
    defineTypeNameAndDebug(binaryWideBand, 0);

    addToRunTimeSelectionTable
    (
        absorptionEmissionModel,
        binaryWideBand,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiationModels::absorptionEmissionModels::binaryWideBand::binaryWideBand
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    absorptionEmissionModel(dict, mesh),
    coeffsDict_(dict.optionalSubDict(typeName + "Coeffs")),
    model1_
    (
        absorptionEmissionModel::New(coeffsDict_.subDict("model1"), mesh)
    ),
    model2_
    (
        absorptionEmissionModel::New(coeffsDict_.subDict("model2"), mesh)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiationModels::absorptionEmissionModels::binaryWideBand::~binaryWideBand()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiationModels::absorptionEmissionModels::binaryWideBand::aCont
(
    const label bandI
) const
{
    return model1_->aCont(bandI) + model2_->aCont(bandI);
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::absorptionEmissionModels::binaryWideBand::aDisp
(
    const label bandI
) const
{
    return model1_->aDisp(bandI) + model2_->aDisp(bandI);
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::absorptionEmissionModels::binaryWideBand::eCont
(
    const label bandI
) const
{
    return model1_->eCont(bandI) + model2_->eCont(bandI);
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::absorptionEmissionModels::binaryWideBand::eDisp
(
    const label bandI
) const
{
    return model1_->eDisp(bandI) + model2_->eDisp(bandI);
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::absorptionEmissionModels::binaryWideBand::ECont
(
    const label bandI
) const
{
    return model1_->ECont(bandI) + model2_->ECont(bandI);
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::absorptionEmissionModels::binaryWideBand::EDisp
(
    const label bandI
) const
{
    return model1_->EDisp(bandI) + model2_->EDisp(bandI);
}

void Foam::radiationModels::absorptionEmissionModels::binaryWideBand::correct
(
    volScalarField& a,
    PtrList<volScalarField>& aLambda
) const
{
    a = dimensionedScalar(dimless/dimLength, 0);
    label nBands = max(model1_->nBands(), model2_->nBands());

    for (label j=0; j<nBands; j++)
    {
        aLambda[j].primitiveFieldRef() = model1_->a(j) + model2_->a(j);

        a.primitiveFieldRef() +=
            aLambda[j].primitiveField();
    }

}
      


// ************************************************************************* //
