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

#ifndef canopyKEpsilonModel_H
#define canopyKEpsilonModel_H

#include "canopyMomentumTransferModel.H"
#include "dynamicMomentumTransportModels.H"
#include "cellSetHashTables.H"

namespace Foam
{

template<class BaseCanopyModel, class CanopykEpsSourceType>
class canopyKEpsilonModel
    : public canopyMomentumTransferModel<BaseCanopyModel>
{

    //- Private Data

      // canopyKEpsilon Dictionary
      dictionary  canopyKEpsilonDict_;

      // TKE Source
      dimensionedScalarCellSet Fk_;

      // Dissipation Source
      dimensionedScalarCellSet Feps_;

      // CanopykEpsSourceModel
      CanopykEpsSourceType kEpsSourceModel_;

      // Reference to k and epsilon
      volScalarField& k_;
      volScalarField& epsilon_;
    
public:
    
    canopyKEpsilonModel
    (
      const treeModel& tree
    );

    ~canopyKEpsilonModel(){};

    virtual wordList addSupFields() const;
    virtual void correctMomentumTransfer();
    virtual dimensionedScalarCellSet& Fturb(const word& name);
    virtual const dimensionedScalarCellSet& Fturb(const word& name) const;

};

}
#ifdef NoRepository
  #include "canopyKEpsilonModel.C"
#endif


#endif
