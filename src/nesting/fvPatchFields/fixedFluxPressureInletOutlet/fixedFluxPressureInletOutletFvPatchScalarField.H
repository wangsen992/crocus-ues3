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
    Foam::fixedFluxPressureInletOutletFvPatchScalarField

Description
    This boundary condition sets the pressure gradient to the provided value
    such that the flux on the boundary is that specified by the velocity
    boundary condition. The pressure gradient is computed in-solver in the same    way fixedFluxPressure works. For outlet, when zeroGradient is applied for 
    velocity, prghHydrostaticPressure is applied using the pressure init. 

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            fixedFluxPressureInletOutlet;
    }
    \endverbatim

See also
    Foam::fixedGradientFvPatchField

SourceFiles
    fixedFluxPressureInletOutletFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef fixedFluxPressureInletOutletFvPatchScalarField_H
#define fixedFluxPressureInletOutletFvPatchScalarField_H

#include "fvPatchFields.H"
#include "mixedFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class fixedFluxPressureInletOutletFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class fixedFluxPressureInletOutletFvPatchScalarField
:
    public mixedFvPatchField<scalar>
{
    // Private Data

        //- Current time index (used for updating)
        label curTimeIndex_;

         //- Name of the velocity field
        word UName_;

        //- Name of the flux transporting the field
        word phiName_;

        //- Name of density field
        word rhoName_;

        //- Name of hydrostatic pressure field
        word ph_rghName_;



public:

    //- Runtime type information
    TypeName("fixedFluxPressureInletOutlet");


    // Constructors

        //- Construct from patch and internal field
        fixedFluxPressureInletOutletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedFluxPressureInletOutletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fixedFluxPressureInletOutletFvPatchScalarField onto
        //  a new patch
        fixedFluxPressureInletOutletFvPatchScalarField
        (
            const fixedFluxPressureInletOutletFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        fixedFluxPressureInletOutletFvPatchScalarField
        (
            const fixedFluxPressureInletOutletFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        fixedFluxPressureInletOutletFvPatchScalarField
        (
            const fixedFluxPressureInletOutletFvPatchScalarField&,
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
                new fixedFluxPressureInletOutletFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        //- Update the patch pressure gradient field from the given snGradp
        virtual void updateCoeffs(const scalarField& snGradp);

        //- Update the patch pressure gradient field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #include "volFields.H"
// 
// namespace Foam
// {
//     template<class GradBC>
//     inline void setSnGrad
//     (
//         volScalarField::Boundary& bf,
//         const FieldField<fvsPatchField, scalar>& snGrad
//     )
//     {
//         forAll(bf, patchi)
//         {
//             if (isA<GradBC>(bf[patchi]))
//             {
//                 refCast<GradBC>(bf[patchi]).updateCoeffs(snGrad[patchi]);
//             }
//         }
//     }
// 
//     template<class GradBC>
//     inline void setSnGrad
//     (
//         volScalarField::Boundary& bf,
//         const tmp<FieldField<fvsPatchField, scalar>>& tsnGrad
//     )
//     {
//         setSnGrad<GradBC>(bf, tsnGrad());
//     }
// }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
