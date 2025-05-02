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
Class
    Foam::canopyAbsorptionEmission

Description
    Provide an absorptionEmissionModel for canopies in the model, by loading
    aggregate radiative properties through objectRegistry

    Provides setting options for wideBand type specification.

    Note: Band 0 must be shortwave and 1 must be longwave
    
SourceFiles
    canopyAbsorptionEmission.C
\*---------------------------------------------------------------------------*/
#ifndef canopyAbsorptionEmission_H
#define canopyAbsorptionEmission_H

#include "absorptionEmissionModel.H"

namespace Foam
{
namespace radiationModels
{
namespace absorptionEmissionModels
{
class canopy
  :
    public absorptionEmissionModel
{
    // Private Data

      //- Coefficients dictionary
      dictionary coeffsDict_;

      //- Reference to the tree-related coefficients
      FixedList<Vector2D<scalar>, 2> iBands_;

      //- Bands
      label nBands_;


public:
    
    //- Runtime type information
    TypeName("canopy");

    // Constructors

      //- Construct from components
      canopy(const dictionary& dict, const fvMesh& mesh);

    //- Destructor
    virtual ~canopy();

    // Member Functions

      //- Absorption coefficient for the dispersed phase
      virtual tmp<volScalarField> aDisp(const label bandI = 0) const;

      //- Emission coefficient for the dispersed phase
      virtual tmp<volScalarField> eDisp(const label bandI = 0) const;

      //- Return emission contribution for the dispersed phase
      virtual tmp<volScalarField> EDisp(const label bandI = 0) const;

      
};

}
}
}

#endif
