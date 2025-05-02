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
\*---------------------------------------------------------------------------*/

#include "canopyKEpsilonModel.H"


namespace Foam
{
template<class BaseCanopyModel, class CanopykEpsSourceType>
canopyKEpsilonModel<BaseCanopyModel, CanopykEpsSourceType>::canopyKEpsilonModel
  (
      const treeModel& tree
  ):
  canopyMomentumTransferModel<BaseCanopyModel>(tree),
  Fk_(BaseCanopyModel::canopyCells().size()),
  Feps_(BaseCanopyModel::canopyCells().size()),
  kEpsSourceModel_
  (
    canopyMomentumTransferModel<BaseCanopyModel>::momentumTransferDict().subDict
    (
      CanopykEpsSourceType::typeName() + "Coeffs"
    )
  ),
  k_(tree.mesh().lookupObjectRef<volScalarField>("k")),
  epsilon_(tree.mesh().lookupObjectRef<volScalarField>("epsilon"))
{
    // Init source cells
    for
    (auto iter = BaseCanopyModel::canopyCells().cbegin(); 
     iter != BaseCanopyModel::canopyCells().cend();
     iter++
    )
    {
        Fk_.set(iter.key(), dimensionedScalar(dimVelocity*dimVelocity/dimTime, 0));
        Feps_.set(iter.key(), dimensionedScalar(dimVelocity*dimVelocity/dimTime/dimTime, 0));
    }

    // Correct momentum transfer
    correctMomentumTransfer();
  
}

template<class BaseCanopyModel, class CanopykEpsSourceType>
wordList canopyKEpsilonModel<BaseCanopyModel, CanopykEpsSourceType>::addSupFields() const
{
    wordList supFields = canopyMomentumTransferModel<BaseCanopyModel>::addSupFields();
    supFields.append(wordList({"k", "epsilon"}));
    return supFields;
}

template<class BaseCanopyModel, class CanopykEpsSourceType>
void canopyKEpsilonModel<BaseCanopyModel, CanopykEpsSourceType>::correctMomentumTransfer()
{
    canopyMomentumTransferModel<BaseCanopyModel>::correctMomentumTransfer();

    // Add turb corrections
    for
    (auto iter = BaseCanopyModel::canopyCells().cbegin(); 
     iter != BaseCanopyModel::canopyCells().cend();
     iter++
    )
    {
        Fk_[iter.key()] = kEpsSourceModel_.Fk
          (
            dimensionedVector
            (
              dimVelocity,
              canopyKEpsilonModel<BaseCanopyModel, CanopykEpsSourceType>::U()[iter.key()]
            ),
            canopyKEpsilonModel<BaseCanopyModel, CanopykEpsSourceType>::Fu()[iter.key()],
            canopyKEpsilonModel<BaseCanopyModel, CanopykEpsSourceType>::lad()[iter.key()],
            canopyKEpsilonModel<BaseCanopyModel, CanopykEpsSourceType>::Cd()
          );
        Feps_[iter.key()] = kEpsSourceModel_.Feps
          (
            dimensionedScalar
            (
              dimVelocity * dimVelocity,
              k_[iter.key()]
            ),
            dimensionedScalar
            (
              dimVelocity * dimVelocity / dimTime,
              epsilon_[iter.key()]
            ),
            dimensionedVector
            (
              dimVelocity,
              canopyKEpsilonModel<BaseCanopyModel, CanopykEpsSourceType>::U()[iter.key()]
            ),
            canopyKEpsilonModel<BaseCanopyModel, CanopykEpsSourceType>::Fu()[iter.key()],
            canopyKEpsilonModel<BaseCanopyModel, CanopykEpsSourceType>::lad()[iter.key()],
            canopyKEpsilonModel<BaseCanopyModel, CanopykEpsSourceType>::Cd()
          );
      }
}

template<class BaseCanopyModel, class CanopykEpsSourceType>
dimensionedScalarCellSet& canopyKEpsilonModel<BaseCanopyModel, CanopykEpsSourceType>::Fturb(const word& name)
{
    if (name == "k")
    {
        return Fk_;
    }
    else if (name == "epsilon")
    {
        return Feps_;
    }
    else
    {
        FatalErrorInFunction 
          << "Input variable " << name << " is not available." 
          << "Available variable names are : k, epsilon" 
          << endl;
    }
                      
}

template<class BaseCanopyModel, class CanopykEpsSourceType>
const dimensionedScalarCellSet& canopyKEpsilonModel<BaseCanopyModel, CanopykEpsSourceType>::Fturb(const word& name) const
{
    if (name == "k")
    {
        return Fk_;
    }
    else if (name == "epsilon")
    {
        return Feps_;
    }
    else
    {
        FatalErrorInFunction 
          << "Input variable " << name << " is not available." 
          << "Available variable names are : k, epsilon" 
          << endl;
    }
                      
}
}
