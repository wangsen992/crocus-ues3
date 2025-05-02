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
    Foam::canopyModel

Description
    Base class for the handling the physics & evapotranspiration of leaves & canopies, 
    which provides source terms to the treeModel that ultimately interacts with 
    the solver. Only to be called by the treeModel. 

    The key is on how leaves spatial domain are defined. With OpenFOAM utilities,
    cellSet can be used (or should ultimately be used), but the given the 
    potential large number of instances, there could be many surface files, or 
    many dictionaries. So this selection should provide multiple methods and 
    potentially easier interface to set up the vegetation coverage. 
    
SourceFiles
    canopyModel.C
\*---------------------------------------------------------------------------*/

#include "canopyModel.H"
// #include "canopySurfaceModel.H"
#include "canopyTriSurfaceModel.H"
#include "canopyKEpsilonModel.H"

#include "canopyEnergyTransferModel.H"

#include "canopyMomentumTransferModel.H"
#include "compressibleMomentumTransportModel.H"
#include "addToRunTimeSelectionTable.H"

#include "HiraokakEpsSourceModel.H"
#include "canopyCellSetModel.H"

namespace Foam
{
// Basic drag model without turbulence source terms
typedef canopyEnergyTransferModel
  <
    canopyMomentumTransferModel
    <
      canopyTriSurfaceModel
      <
        canopyModel
      >    
    >
  > 
    triSurfaceDragCanopyModel;

addNamedToRunTimeSelectionTable
(
  canopyModel, 
  triSurfaceDragCanopyModel, 
  treeModel, 
  triSurfaceDragCanopyModel
);

// Drag model with Hiraoka model for RAS model correction
typedef
  canopyEnergyTransferModel
  <
    canopyKEpsilonModel
    <
      canopyTriSurfaceModel
      <
        canopyModel
      >, 
      HiraokakEpsSourceModel
    >
  > 
    triSurfaceKEpsilonCanopyModel;

addNamedToRunTimeSelectionTable
(
  canopyModel, 
  triSurfaceKEpsilonCanopyModel, 
  treeModel, 
  triSurfaceKEpsilonCanopyModel
);

typedef canopyEnergyTransferModel
  <
    canopyMomentumTransferModel
    <
      canopyCellSetModel
      <
        canopyModel
      >    
    >
  > 
    cellSetDragCanopyModel;

addNamedToRunTimeSelectionTable
(
  canopyModel, 
  cellSetDragCanopyModel, 
  treeModel, 
  cellSetDragCanopyModel
);

typedef
  canopyEnergyTransferModel
  <
    canopyKEpsilonModel
    <
      canopyCellSetModel
      <
        canopyModel
      >, 
      HiraokakEpsSourceModel
    >
  > 
    cellSetKEpsilonCanopyModel;

addNamedToRunTimeSelectionTable
(
  canopyModel, 
  cellSetKEpsilonCanopyModel,
  treeModel, 
  cellSetKEpsilonCanopyModel
);
}
