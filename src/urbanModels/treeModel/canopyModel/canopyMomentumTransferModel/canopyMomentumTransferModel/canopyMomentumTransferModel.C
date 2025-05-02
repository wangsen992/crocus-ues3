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

template<class BaseCanopyModel>
void Foam::canopyMomentumTransferModel<BaseCanopyModel>::readCoeffs()
{
    Cd_ = momentumTransferDict_.lookup<double>("Cd");
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
    Fu_(BaseCanopyModel::canopyCells().size())
{
    Info << "Init in canopyMomentumTransferModel.." << endl;
    readCoeffs();
    correctMomentumTransfer();
};

template<class BaseCanopyModel>
void Foam::canopyMomentumTransferModel<BaseCanopyModel>::correctMomentumTransfer(){
    vector Ui;
    scalar magUi;
    vector ladi;
    label celli;

    const labelHashSet& cells = BaseCanopyModel::canopyCells();
    forAll(cells, i)
    {
        celli = cells.sortedToc()[i];
        Ui = U_[celli];
        magUi = mag(Ui);
        ladi = BaseCanopyModel::lad()[celli].value();

        // Update drag force term
        Fu_.set
        (celli, 
         dimensionedVector
         (
            dimVelocity/dimTime,
            - Cd_ * cmptMultiply(ladi, magUi * Ui)
         )
        );
    }
}


template<class BaseCanopyModel>
dimensionedScalarCellSet& Foam::canopyMomentumTransferModel<BaseCanopyModel>::Fturb(const word& name)
{
    NotImplemented;
}
template<class BaseCanopyModel>
const dimensionedScalarCellSet& Foam::canopyMomentumTransferModel<BaseCanopyModel>::Fturb(const word& name) const
{
    NotImplemented;
}
