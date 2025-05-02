/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "fixedFluxPressureInletOutletFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fixedFluxPressureInletOutletFvPatchScalarField::fixedFluxPressureInletOutletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(p, iF),
    curTimeIndex_(-1),
    UName_("U"),
    phiName_("phi"),
    rhoName_("rho"),
    ph_rghName_("ph_rgh")
{}


Foam::fixedFluxPressureInletOutletFvPatchScalarField::fixedFluxPressureInletOutletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<scalar>(p, iF, dict),
    curTimeIndex_(-1),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    ph_rghName_(dict.lookupOrDefault<word>("ph_rgh", "ph_rgh"))
{
    if (dict.found("value") && dict.found("gradient"))
    {
        fvPatchField<scalar>::operator=(scalarField("value", dict, p.size()));
        this->refGrad() = scalarField("gradient", dict, p.size());
    }
    else
    {
        fvPatchField<scalar>::operator=(patchInternalField());
        this->refGrad() = Zero;
    }
}


Foam::fixedFluxPressureInletOutletFvPatchScalarField::fixedFluxPressureInletOutletFvPatchScalarField
(
    const fixedFluxPressureInletOutletFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<scalar>(ptf, p, iF, mapper),
    curTimeIndex_(-1),
    UName_(ptf.UName_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    ph_rghName_(ptf.ph_rghName_)
{}


Foam::fixedFluxPressureInletOutletFvPatchScalarField::fixedFluxPressureInletOutletFvPatchScalarField
(
    const fixedFluxPressureInletOutletFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(ptf, iF),
    curTimeIndex_(-1),
    UName_(ptf.UName_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    ph_rghName_(ptf.ph_rghName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fixedFluxPressureInletOutletFvPatchScalarField::updateCoeffs
(
    const scalarField& snGradp
)
{
    if (updated())
    {
        return;
    }
    
    // Update 
      const scalarField& rhop =
        patch().lookupPatchField<volScalarField, scalar>(rhoName_);

    const scalarField& ph_rghp =
        patch().lookupPatchField<volScalarField, scalar>(ph_rghName_);

    const scalarField& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);

    const vectorField& Up =
        patch().lookupPatchField<volVectorField, vector>(UName_);

    this->refValue()=
    (
        ph_rghp
      - 0.5*rhop*(1.0 - pos0(phip))*magSqr(Up)
    );

    vector Utarget(3, 3, 0);
    scalar t = this->db().time().value();
	  // if (t > 100)
	  // {
	  // 	// Utarget = vector(5 * sin(6.28*t/1000) + 1 * sin(6.28*t/300), 5 * cos(6.28*t/2000) + 0.5 * cos(6.28*t/450),0);
	  // 	Utarget = 100/100 * vector(5 * sin(6.28*t/1000) + 1 * sin(6.28*t/3000), 5 * cos(6.28*t/2000) + 0.5 * cos(6.28*t/4500),0);
	  // }
	  // else 
	  // {
	  // 	// Utarget = 100/100 * vector(5 * sin(6.28*t/1000) + 1 * sin(6.28*t/300), 5 * cos(6.28*t/2000) + 0.5 * cos(6.28*t/450),0);
	  // 	Utarget = t/100 * vector(5 * sin(6.28*t/1000) + 1 * sin(6.28*t/3000), 5 * cos(6.28*t/2000) + 0.5 * cos(6.28*t/4500),0);
	  // }

    // this->valueFraction() = pos0(phip);
    this->valueFraction() = pos0(Utarget & this->patch().Sf());

    // Update index and gradient
    curTimeIndex_ = this->db().time().timeIndex();

    this->refGrad() = snGradp;
    mixedFvPatchField<scalar>::updateCoeffs();
}


void Foam::fixedFluxPressureInletOutletFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }


    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        FatalErrorInFunction
            << "updateCoeffs(const scalarField& snGradp) MUST be called before"
               " updateCoeffs() or evaluate() to set the boundary gradient."
            << exit(FatalError);
    }
}


void Foam::fixedFluxPressureInletOutletFvPatchScalarField::write(Ostream& os) const
{
    mixedFvPatchField<scalar>::write(os);
    writeEntryIfDifferent<word>(os, "U", "U", UName_);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);
    writeEntryIfDifferent<word>(os, "ph_rgh", "ph_rgh", ph_rghName_);
    // writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        fixedFluxPressureInletOutletFvPatchScalarField
    );
}


// ************************************************************************* //
