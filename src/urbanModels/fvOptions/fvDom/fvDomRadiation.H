/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2021 OpenFOAM Foundation
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
    Foam::fv::fvDomRadiation

Description
    Calculates and applies the fvDomRadiation source to the energy equation.

Usage
    Example usage:
    \verbatim
    radiationCoeffs
    {
        libs    ("libradiationModels.so");
        type    fvDomRadiation;
    }
    \endverbatim

SourceFiles
    fvDomRadiation.C

\*---------------------------------------------------------------------------*/

#ifndef fvDomRadiation_H
#define fvDomRadiation_H

#include "fvModel.H"
#include "uniformDimensionedFields.H"
#include "radiationModel.H"
#include "fvDOM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                          Class fvDomRadiation Declaration
\*---------------------------------------------------------------------------*/

class fvDomRadiation
:
    public fvModel
{
    // Private Data

        //- The fvDomRadiation model pointer
        mutable autoPtr<radiationModels::fvDOM> radiation_;


public:

    //- Runtime type information
    TypeName("fvDomRadiation");


    // Constructors

        //- Construct from explicit source name and mesh
        fvDomRadiation
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Disallow default bitwise copy construction
        fvDomRadiation(const fvDomRadiation&) = delete;


    // Member Functions

        // Checks

            //- Return the list of fields for which the option adds source term
            //  to the transport equation
            virtual wordList addSupFields() const;


        // Evaluate

            //- Add explicit contribution to compressible energy equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<scalar>& eqn,
                const word& fieldName
            ) const;

            //- Add explicit contribution to compressible energy equation
            virtual void addSup
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                fvMatrix<scalar>& eqn,
                const word& fieldName
            ) const;

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const fvDomRadiation&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
