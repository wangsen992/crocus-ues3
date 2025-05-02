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
    Foam::canopyMomentumTransferModel

Description
    Base class to provide momentumTransport (include velocity, and associated turbulence
    parameters. 
    
SourceFiles
    canopyMomentumTransferModel.C
\*---------------------------------------------------------------------------*/

#ifndef canopyMomentumTransferModel_H
#define canopyMomentumTransferModel_H

#include "treeModel.H"
#include "dynamicMomentumTransportModel.H"

// typedefs

typedef compressible::momentumTransportModel momentumTransportModel;

namespace Foam
{

template<class BaseCanopyModel>
class canopyMomentumTransferModel
:
    public BaseCanopyModel
{
    
    // private data
protected:

        //- constant drag coffient
        scalar Cd_;

        //- Dictionary for momentumTransfer
        dictionary momentumTransferDict_;

        //- Reference to velocity
        volVectorField& U_;

        //- Holding values for the drag force term
        dimensionedVectorCellSet Fu_;


    // Private Function
    
        void readCoeffs();

public:

    canopyMomentumTransferModel
    (
      const treeModel& tree
    );

    virtual ~canopyMomentumTransferModel(){};

    //- Return a list of field names that will be impacted
    virtual wordList addSupFields() const {return {"U"};};

    const volVectorField& U() const {return U_;}
    scalar Cd(){return Cd_;}
    dictionary& momentumTransferDict(){return momentumTransferDict_;}
    virtual dimensionedVectorCellSet& Fu(){return Fu_;};
    virtual dimensionedScalarCellSet& Fturb(const word&);
    virtual const dimensionedScalarCellSet& Fturb(const word&) const;
    virtual const dimensionedVectorCellSet& Fu() const {return Fu_;};
    virtual void correctMomentumTransfer();

    
};

}

#ifdef NoRepository
  #include "canopyMomentumTransferModel.C"
#endif

#endif

