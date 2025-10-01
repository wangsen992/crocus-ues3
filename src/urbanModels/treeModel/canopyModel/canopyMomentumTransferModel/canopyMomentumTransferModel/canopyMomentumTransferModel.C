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
\*---------------------------------------------------------------------------*/

#include "canopyMomentumTransferModel.H"
#include "dimensionSets.H"

template<class BaseCanopyModel>
void Foam::canopyMomentumTransferModel<BaseCanopyModel>::readCoeffs()
{
    Cd_ = momentumTransferDict_.get<scalar>("Cd");
}

template<class BaseCanopyModel>
Foam::canopyMomentumTransferModel<BaseCanopyModel>::canopyMomentumTransferModel
(
    const treeModel& tree
)
:
    BaseCanopyModel(tree),
    momentumTransferDict_(BaseCanopyModel::dict().subDict("momentumTransfer")),
    Cd_(),
    U_(tree.mesh().lookupObjectRef<volVectorField>("U")),
    Fu_
    (
      IOobject
      (
        "canopy_Fu",
        tree.mesh().time().constant(),
        tree.mesh(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
      ),
      tree.mesh(),
      dimensionedVector(dimVelocity/dimTime, vector(0,0,0))
    )
{
    Info << "Init in canopyMomentumTransferModel.." << endl;
    readCoeffs();
    correctMomentumTransfer();
};

template<class BaseCanopyModel>
void Foam::canopyMomentumTransferModel<BaseCanopyModel>::correctMomentumTransfer(){
        // Update drag force term
        Fu_ = - Cd_ * BaseCanopyModel::lad() * mag(U_) * U_;
}


template<class BaseCanopyModel>
volScalarField& Foam::canopyMomentumTransferModel<BaseCanopyModel>::Fturb(const word& name)
{
    NotImplemented;
}
template<class BaseCanopyModel>
const volScalarField& Foam::canopyMomentumTransferModel<BaseCanopyModel>::Fturb(const word& name) const
{
    NotImplemented;
}
