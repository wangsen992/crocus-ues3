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

#include "solarRadiationMixedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

#include "fvDOM.H"
#include "wideBand.H"
#include "constants.H"

using namespace Foam::constant;
using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solarRadiationMixedFvPatchScalarField::
solarRadiationMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    radiationCoupledBase(p, "undefined", scalarField::null()),
    TName_("T"),
    thetaStart_(),
    thetaRange_(),
    deltaTheta_(),
    phi0_()
{
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 1.0;
}


Foam::solarRadiationMixedFvPatchScalarField::
solarRadiationMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    radiationCoupledBase(p, dict),
    TName_(dict.lookupOrDefault<word>("T", "T")),
    thetaStart_(dict.lookup<scalar>("thetaStart")),
    thetaRange_(dict.lookup<scalar>("thetaRange")),
    deltaTheta_(dict.lookup<scalar>("deltaTheta")),
    phi0_(dict.lookup<scalar>("phi0"))
{
    if (dict.found("value"))
    {
        fvPatchScalarField::operator=
        (
            scalarField("value", dict, p.size())
        );
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        refValue() = 0;
        refGrad() = 0.0;

        fvPatchScalarField::operator=(refValue());
    }
}


Foam::solarRadiationMixedFvPatchScalarField::
solarRadiationMixedFvPatchScalarField
(
    const solarRadiationMixedFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    radiationCoupledBase
    (
        p,
        ptf.emissivityMethod(),
        ptf.emissivity_,
        mapper
    ),
    TName_(ptf.TName_),
    thetaStart_(ptf.thetaStart_),
    thetaRange_(ptf.thetaRange_),
    deltaTheta_(ptf.deltaTheta_),
    phi0_(ptf.phi0_)
{}


Foam::solarRadiationMixedFvPatchScalarField::
solarRadiationMixedFvPatchScalarField
(
    const solarRadiationMixedFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    radiationCoupledBase
    (
        ptf.patch(),
        ptf.emissivityMethod(),
        ptf.emissivity_
    ),
    TName_(ptf.TName_),
    thetaStart_(ptf.thetaStart_),
    thetaRange_(ptf.thetaRange_),
    deltaTheta_(ptf.deltaTheta_),
    phi0_(ptf.phi0_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::solarRadiationMixedFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchScalarField::autoMap(m);
    radiationCoupledBase::autoMap(m);
}


void Foam::solarRadiationMixedFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);
    radiationCoupledBase::rmap(ptf, addr);
}


void Foam::solarRadiationMixedFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;

    const radiationModels::fvDOM& dom =
        db().lookupObject<radiationModels::fvDOM>("radiationProperties");

    label rayId = -1;
    label lambdaId = -1;
    dom.setRayIdLambdaId(internalField().name(), rayId, lambdaId);

    // Load dom aggregate info
    vectorField dList(dom.nRay());
    vectorField dAveList(dom.nRay());
    scalarField omegaList(dom.nRay());
    DynamicList<label> downRayList(dom.nRay());
    DynamicList<label> upRayList(dom.nRay());

    // Separate rays into upward and downward
    for (label i = 0; i < dom.nRay(); i++)
    {
      dList[i] = dom.IRay(i).d();
      dAveList[i] = dom.IRay(i).dAve();
      omegaList[i] = dom.IRay(i).omega();
      if ((dom.IRay(i).dAve() & vector(0,0,1)) > 0.0)
      {
          upRayList.append(i);
      }
      else
      {
          downRayList.append(i);
      }
    }
    upRayList.shrink();
    downRayList.shrink();

    // Use solar direction calc
    // Set 
    scalar theta_start = thetaStart_;

    // For a normal day simulation, deltaTheta = pi / 43200

    scalar theta0 = theta_start - deltaTheta_ * db().time().timeOutputValue();
    vector d0 = vector(-sin(theta0)*cos(phi0_), sin(phi0_), -cos(theta0));
    if (mag(d0.z()) < VSMALL)
    {
        d0.z() += sign(d0.z()) * VSMALL;
    }

    // Find aligned angle to direct beam
    scalarField dprodList(dList & d0);
    label i0 = findMax(dprodList);
    vector dAve0 = dAveList[i0];
    scalar nAve0 = mag(dAve0 & vector(0,0,1));
    scalar nAveDown = sum(mag(vectorField(dAveList, downRayList) & vector(0,0,1)));

    Info << "cos(theta0) = " << cos(theta0) << endl;
    
    // Sky radiation properties
    scalar I0_0 = 0.0;
    scalar I0_1 = 0.0;
    scalar Idefault_0 = 0.0;
    scalar Idefault_1 = 0.0;
    if ( d0.z() < 0) // daytime condition
    {
      const scalar R_DV = 600 * exp(-0.185 / cos(theta0)) * cos(theta0);
      const scalar R_dV = 0.4 * (600 - R_DV) * cos(theta0);
      const scalar w = 1320 * pow(10, -1.1950 + 0.4459 * log10(1/cos(theta0)) - 0.0345 * sqr(log10(1/cos(theta0))));
      const scalar R_DN = (720 * exp(-0.06 / cos(theta0)) - w) * cos(theta0);
      const scalar R_dN = 0.6 * (720 - R_DN - w) * cos(theta0);

      I0_0 = R_DV / nAve0;
      I0_1 = R_DN / nAve0;
      Idefault_0 = R_dV / nAveDown;
      Idefault_1 = R_dN / nAveDown;
    }
    else // night-time condition
    {
      // Temporary setting for night radiation
      const scalar R_dN = 150;
      // Night time diffusive longwave to reduce the cooling speed
      Idefault_1 = R_dN / nAveDown;
    }

    // End of solar calc
    const label patchi = patch().index();

    if (dom.nLambda() == 0)
    {
        FatalErrorInFunction
            << " a non-grey boundary condition is used with a grey "
            << "absorption model" << nl << exit(FatalError);
    }

    scalarField& Iw = *this;
    const vectorField n(patch().Sf()/patch().magSf());

    radiationModels::radiativeIntensityRay& ray =
        const_cast<radiationModels::radiativeIntensityRay&>(dom.IRay(rayId));

    const scalarField nAve(n & ray.dAve());
    Info << "[Debug] dAve= " << ray.dAve() << endl;

    ray.qr().boundaryFieldRef()[patchi] += Iw*nAve;

    const scalarField Eb
    (
        dom.blackBody().bLambda(lambdaId).boundaryField()[patchi]
    );

    scalarField temissivity = emissivity();

    scalarField& qem = ray.qem().boundaryFieldRef()[patchi];
    scalarField& qin = ray.qin().boundaryFieldRef()[patchi];

    // Use updated Ir while iterating over rays
    // avoids to used lagged qin
    scalarField Ir = dom.IRay(0).qin().boundaryField()[patchi];

    for (label rayI=1; rayI < dom.nRay(); rayI++)
    {
        Ir += dom.IRay(rayI).qin().boundaryField()[patchi];
    }

    if (lambdaId == 0)
    {
        if (rayId == i0)
        {
            Info << "Update direct: " << I0_0 
                 << ", " << max(Iw) << endl;
            refValue() = I0_0;
            valueFraction() = 1.0;
            refGrad() = 0.0;
            qem += refValue() * nAve;
        }
        else if (findIndex(downRayList, rayId) != -1)
        {
            Info << "Update diffusive: " << Idefault_0 
                 << ", " << max(Iw) << endl;
            refValue() = Idefault_0;
            valueFraction() = 1.0;
            refGrad() = 0.0;
            qem += refValue() * nAve;
        }
        else
        {
            valueFraction() = 0.0;
            refGrad() = 0.0;
            refValue() = 0.0; // not used
            qin = 0;
        }
    }
    else if (lambdaId == 1)
    {
        if (rayId == i0)
        {
            // Treat the main direction as the else condition here to avoid
            // wrongly assigned value
            if (d0.z() < 0)
            {
              refValue() = I0_1;
              valueFraction() = 1.0;
              refGrad() = 0.0;
              qem += refValue() * nAve;
            }
            else
            {
              refValue() = 0.0;
              valueFraction() = 0.0;
              refGrad() = 0.0;
              qin = 0;
            }
        }
        else if (findIndex(downRayList, rayId) != -1)
        {
            refValue() = Idefault_1;
            valueFraction() = 1.0;
            refGrad() = 0.0;
            qem += refValue() * nAve;
        }
        else
        {
            refValue() = 0.0;
            valueFraction() = 0.0;
            refGrad() = 0.0;
            qin = 0;
        }
    }


    // if (this->patch().name() == "floor")
    // {
    // Info << "average(Iw) = " << average(Iw) << ", " << rayId << ", " << lambdaId << ", " << dom.IRay(rayId).d() << endl;
    // }
    // Restore tag
    UPstream::msgType() = oldTag;

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::solarRadiationMixedFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    radiationCoupledBase::write(os);
    writeEntryIfDifferent<word>(os, "T", "T", TName_);
    writeEntry(os, "thetaStart", thetaStart_);
    writeEntry(os, "thetaRange", thetaRange_);
    writeEntry(os, "deltaTheta", deltaTheta_);
    writeEntry(os, "phi0", phi0_);

    // [TODO] include writeEntry for the required parameters to restart case
    // e.g. thetaStart, thetaRange, deltaTheta, etc. 
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        solarRadiationMixedFvPatchScalarField
    );
}


// ************************************************************************* //
