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
    Foam::canopyKEpsilonModel

Description
    Base class to provide momentumTransport (include velocity, and associated turbulence
    parameters. 
    
SourceFiles
    canopyKEpsilonModel.C
\*---------------------------------------------------------------------------*/
#ifndef HiraokakEpsSourceModel_H
#define HiraokakEpsSourceModel_H

#include "fvCFD.H"

namespace Foam
{

class HiraokakEpsSourceModel
{
    
    // Private Data
        
        // Coefficients
        scalar  C_peps1_;

public: 

    HiraokakEpsSourceModel
    (
      const dictionary& dict
    );
    
    static word typeName()
    {
        return "Hiraoka";
    }

    dimensionedScalar Fk
    (
      dimensionedVector U,
      dimensionedVector Fu,
      dimensionedVector lad,
      scalar Cd
    )
    {
        return U & ((-1) * Fu);
    }

    dimensionedScalar Feps
    (
      dimensionedScalar k,
      dimensionedScalar eps,
      dimensionedVector U,
      dimensionedVector Fu,
      dimensionedVector lad,
      scalar Cd
    )
    {
        return eps / k * C_peps1_ * pow(k, 1.5) * mag(lad);
    }
};
}
#endif
