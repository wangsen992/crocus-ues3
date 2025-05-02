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

Class
    Foam::energyBalanceHeatFluxTemperatureFvPatchScalarField

Description
    
    This is a adapted boundary condition, no coupling available. 


    This boundary condition applies a heat flux condition to temperature
    on an external wall in one of three modes:

      - fixed power: supply Q
      - fixed heat flux: supply q
      - fixed heat transfer coefficient: supply h and Ta

    where:
    \vartable
        Q  | Power [W]
        q  | Heat flux [W/m^2]
        h  | Heat transfer coefficient [W/m^2/K]
        Ta | Ambient temperature [K]
    \endvartable

    For heat transfer coefficient mode optional thin thermal layer resistances
    can be specified through thicknessLayers and kappaLayers entries.

    The thermal conductivity \c kappa can either be retrieved from various
    possible sources, as detailed in the class temperatureCoupledBase.

    The ambient temperature Ta is specified as a Foam::Function1 of time but
    uniform is space.

Usage
    \table
    Property     | Description                 | Required | Default value
    qr           | Name of the radiative field | no | none
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            energyBalanceHeatFluxTemperature;

        value           $internalField;
    }
    \endverbatim

See also
    Foam::mixedFvPatchScalarField

SourceFiles
    energyBalanceHeatFluxTemperatureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef energyBalanceHeatFluxTemperatureFvPatchScalarField_H
#define energyBalanceHeatFluxTemperatureFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "temperatureCoupledBase.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class energyBalanceHeatFluxTemperatureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class energyBalanceHeatFluxTemperatureFvPatchScalarField
:
    public mixedFvPatchScalarField
{
public:

    // Public data

        //- Operation mode enumeration
        enum operationMode
        {
            fixedPower,
            fixedHeatFlux,
            fixedHeatTransferCoeff
        };

        static const NamedEnum<operationMode, 3> operationModeNames;


private:

    // Private Data



public:

    //- Runtime type information
    TypeName("energyBalanceHeatFluxTemperature");


    // Constructors

        //- Construct from patch and internal field
        energyBalanceHeatFluxTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        energyBalanceHeatFluxTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  energyBalanceHeatFluxTemperatureFvPatchScalarField
        //  onto a new patch
        energyBalanceHeatFluxTemperatureFvPatchScalarField
        (
            const energyBalanceHeatFluxTemperatureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        energyBalanceHeatFluxTemperatureFvPatchScalarField
        (
            const energyBalanceHeatFluxTemperatureFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        energyBalanceHeatFluxTemperatureFvPatchScalarField
        (
            const energyBalanceHeatFluxTemperatureFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new energyBalanceHeatFluxTemperatureFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        // Access

            //- Allow manipulation of the boundary values
            virtual bool fixesValue() const
            {
                return false;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            //  Used to update fields following mesh topology change
            virtual void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            //  Used to reconstruct fields
            virtual void rmap(const fvPatchScalarField&, const labelList&);


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // I-O

            //- Write
            void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
