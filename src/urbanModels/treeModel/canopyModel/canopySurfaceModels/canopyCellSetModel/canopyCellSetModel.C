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
    Foam::canopyGeometryModel

\*---------------------------------------------------------------------------*/

#include "canopyCellSetModel.H"
#include "dimensionedScalarFwd.H"
#include "triSurfaceSearch.H"
#include "pointIndexHit.H"
#include "surfaceMeshTools.H"
#include "ListOps.H"
#include "cell.H"

#include "IOstreams.H"

namespace Foam
{

template<class BaseCanopyModel>
void canopyCellSetModel<BaseCanopyModel>::init()
{
    const polyMesh& mesh = BaseCanopyModel::tree().mesh();

    // Extract cellSet for canopy cells
    labelHashSet canopyCellsIndex(set_);
    
    if (Pstream::parRun())
    {
      Info << "Parallel Run..." << endl;
      Info << "Pstream::masterNo() =  "  << Pstream::masterNo() <<  endl;
      Pout << "Pstream::myProcNo() =  "  << Pstream::myProcNo() <<  endl;
      Pout << "mesh.cells().size() = "  << mesh.cells().size() << endl;
      Pout <<  "canopyCellsIndex.size() = " << canopyCellsIndex.size() << endl;
    }
    else
    {
      Info << "Serial Run..." << endl;
    }

    this->canopyCells() = canopyCellsIndex;

    // Init properties
    HashTable<vector, label> cellLeafArea(canopyCellsIndex.size());
    HashTable<scalarField, label> cellDiameterList(canopyCellsIndex.size());

    forAll(canopyCellsIndex, celli)
    {
        cellLeafArea.set(canopyCellsIndex.toc()[celli], vector(0,0,0));
        cellDiameterList.set(canopyCellsIndex.toc()[celli], scalarField());
    }

    labelList cellList = cellLeafArea.toc();

    // constant properties assigned to all tree cells
    dimensionedScalar la(dimArea/dimVolume, this->surfaceModelDict().lookupOrDefault("la", 0.2));
    dimensionedScalar lad(dimArea/dimVolume, this->surfaceModelDict().lookupOrDefault("lad", 1));
    dimensionedScalar ldia(dimArea/dimVolume, this->surfaceModelDict().lookupOrDefault("ldia", 1));
    dimensionedScalar laLit(dimArea/dimVolume, this->surfaceModelDict().lookupOrDefault("laLit", 1));
    dimensionedScalar laCov(dimArea/dimVolume, this->surfaceModelDict().lookupOrDefault("laCov",1));
    forAll(cellList, i)
    {
        this->la().set(cellList[i], la);
        this->lad().set(cellList[i], lad);
        this->ldia().set(cellList[i], ldia);
        this->laLit().set(cellList[i], laLit);
        this->laCov().set(cellList[i], laCov);
    }
}

template<class BaseCanopyModel>
canopyCellSetModel<BaseCanopyModel>::canopyCellSetModel
(
    const treeModel& tree
)
:
    canopySurfaceModel<BaseCanopyModel>(tree),
    set_(tree.mesh(), "c0") // load the cellSet
{
    init();
    // Testing find lalit
    vector defaultDirection(0, 0, -1);
};
}

template<class BaseCanopyModel>
dimensionedScalarCellSet canopyCellSetModel<BaseCanopyModel>::correctLaCov
(
  vector direction
)
{
    Info << "Update laCov with no change" << direction << endl;
    return this->laCov();

} 

