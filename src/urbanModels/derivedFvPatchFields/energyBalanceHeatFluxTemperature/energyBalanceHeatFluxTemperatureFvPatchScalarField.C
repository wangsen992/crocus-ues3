/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "energyBalanceHeatFluxTemperatureFvPatchScalarField.H"
#include "mathematicalConstants.H"
#include "volFields.H"
#include "physicoChemicalConstants.H"
#include "addToRunTimeSelectionTable.H"

using Foam::constant::physicoChemical::sigma;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::energyBalanceHeatFluxTemperatureFvPatchScalarField::
energyBalanceHeatFluxTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF)
{
    refValue() = 0;
    refGrad() = 0;
    valueFraction() = 1;
}


Foam::energyBalanceHeatFluxTemperatureFvPatchScalarField::
energyBalanceHeatFluxTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF)
{

    if(dict.found("value"))
    {
      fvPatchScalarField::operator=(scalarField("value", dict, p.size()));
    }
    else
    {
      fvPatchScalarField::operator=(this->patchInternalField());
    }
    
    
    if (dict.found("refValue"))
    {
        // Full restart
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        refValue() = *this;
        refGrad() = 0;
        valueFraction() = 1;
    }
}


Foam::energyBalanceHeatFluxTemperatureFvPatchScalarField::
energyBalanceHeatFluxTemperatureFvPatchScalarField
(
    const energyBalanceHeatFluxTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper)
{
}


Foam::energyBalanceHeatFluxTemperatureFvPatchScalarField::
energyBalanceHeatFluxTemperatureFvPatchScalarField
(
    const energyBalanceHeatFluxTemperatureFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(tppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::energyBalanceHeatFluxTemperatureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchScalarField::autoMap(m);
}


void Foam::energyBalanceHeatFluxTemperatureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);
}


void Foam::energyBalanceHeatFluxTemperatureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Load related models 
    const scalarField& Tp(*this);

    // Store current valueFraction and refValue for relaxation
    const scalarField valueFraction0(valueFraction());
    const scalarField refValue0(refValue());

    scalarField qr(Tp.size(), 0);
    scalarField qin(Tp.size(), 0);
    scalarField qem(Tp.size(), 0);
    scalarField alphat(Tp.size(), 0);
    scalarField rho(Tp.size(), 0);
    scalarField Cp(Tp.size(), 1000);
    scalarField qEff(Tp.size(), 0);

    qr  = patch().lookupPatchField<volScalarField, scalar> ( "qr");
    qin = patch().lookupPatchField<volScalarField, scalar> ( "qin");
    qem = patch().lookupPatchField<volScalarField, scalar> ( "qem");
    alphat = 
        patch().lookupPatchField<volScalarField, scalar>
        (
          IOobject::groupName("alphat", internalField().group())
        );
    // rho = 
    //     patch().lookupPatchField<volScalarField, scalar>
    //     (
    //       IOobject::groupName("rho", internalField().group())
    //     );
    // Cp = 
    //     patch().lookupPatchField<volScalarField, scalar>
    //     (
    //       IOobject::groupName("Cp", internalField().group())
    //     );
    qEff = Cp  * alphat
          * 
            (this->refValue() - this->patchInternalField()) 
          * this->patch().deltaCoeffs();

    // Temperatue change of the patch is computed by having net energy change /
    // specific heat capacity
    // qr: total radiative heat flux [W/m2]
    refGrad() = 0.0;
    refValue() += (- qr - qEff) / (2400 * 0.1 * 880);
    // refValue() += (qin + qem) / (2400 * 0.1 * 880);
    valueFraction() = 1.0;

    mixedFvPatchScalarField::updateCoeffs();

    if (debug)
    {
    }
}


void Foam::energyBalanceHeatFluxTemperatureFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);

    // writeEntry(os, "refValue", refValue());
    // writeEntry(os, "refGradient", refGrad());
    // writeEntry(os, "valueFraction", valueFraction());
    // writeEntry(os, "value", *this);
    // os.writeEntry("value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        energyBalanceHeatFluxTemperatureFvPatchScalarField
    );
}

// ************************************************************************* //
