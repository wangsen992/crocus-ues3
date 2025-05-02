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
    Foam::radiationModels::absorptionEmissionModels::binaryWideBand

Description
    Radiation coefficient based on two absorption models

SourceFiles
    binaryWideBand.C

\*---------------------------------------------------------------------------*/

#ifndef binary_H
#define binary_H

#include "absorptionEmissionModel.H"
#include "Vector2D.H"
#include "Table.H"
#include "dimensionedScalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{
namespace absorptionEmissionModels
{

/*---------------------------------------------------------------------------*\
                           Class binaryWideBand Declaration
\*---------------------------------------------------------------------------*/

class binaryWideBand
:
    public absorptionEmissionModel
{
    // Private Data

        //- Coefficients dictionary
        dictionary coeffsDict_;

        //- First absorption model
        autoPtr<absorptionEmissionModel> model1_;

        //- Second absorption model
        autoPtr<absorptionEmissionModel> model2_;


public:

    //- Runtime type information
    TypeName("binaryWideBand");


    // Constructors

        //- Construct from components
        binaryWideBand(const dictionary& dict, const fvMesh& mesh);


    //- Destructor
    virtual ~binaryWideBand();


    // Member Functions

        // Access

            // Absorption coefficient

                //- Absorption coefficient for continuous phase
                virtual tmp<volScalarField> aCont(const label bandI = 0) const;

                //- Absorption coefficient for dispersed phase
                virtual tmp<volScalarField> aDisp(const label bandI = 0) const;


            // Emission coefficient

                //- Emission coefficient for continuous phase
                virtual tmp<volScalarField> eCont(const label bandI = 0) const;

                //- Emission coefficient for dispersed phase
                virtual tmp<volScalarField> eDisp(const label bandI = 0) const;


            // Emission contribution

                //- Emission contribution for continuous phase
                virtual tmp<volScalarField> ECont(const label bandI = 0) const;

                //- Emission contribution for continuous phase
                virtual tmp<volScalarField> EDisp(const label bandI = 0) const;

        //- Lower and uppwer limit fo band i
        // Note: This must be included so blackbody module can work
        inline const Vector2D<scalar>& bands(const label bandi) const
        {
          if (model1_->bands(bandi) != Vector2D<scalar>::one)
          {
              return model1_->bands(bandi);
          }
          else if (model2_->bands(bandi) != Vector2D<scalar>::one)
          {
              return model2_->bands(bandi);
          }
          else
          {
              return Vector2D<scalar>::one;
          }
        }

        //- Const access to the number of bands 
        virtual label nBands() const 
        {return max(model2_->nBands(), model1_->nBands());}
        
        //- Correct rays
        virtual void correct
        (
            volScalarField& a,
            PtrList<volScalarField>& aLambda
        ) const;
};
//


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace absorptionEmissionModels
} // End namespace radiationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
