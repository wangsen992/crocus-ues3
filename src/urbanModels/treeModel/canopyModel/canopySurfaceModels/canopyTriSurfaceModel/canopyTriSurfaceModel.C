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

#include "fileName.H"
#include "canopyTriSurfaceModel.H"
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
void canopyTriSurfaceModel<BaseCanopyModel>::init()
{
    const polyMesh& mesh = BaseCanopyModel::tree().mesh();
    // Compute the total area of leaves within a mesh cell
    const pointField& surfacePoints = surface_.points();
    const faceList& faces = surface_.faces();
    // Correct the intersected cells by the leafy surface
    labelHashSet canopyCellsIndex = surfaceMeshTools::findSurfaceCutCells(mesh, surface_);
    
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
    HashTable<vector, label> cellLeafArea(canopyCellsIndex.size());
    HashTable<scalarField, label> cellDiameterList(canopyCellsIndex.size());

    forAll(canopyCellsIndex, celli)
    {
        cellLeafArea.set(canopyCellsIndex.toc()[celli], vector(0,0,0));
        cellDiameterList.set(canopyCellsIndex.toc()[celli], scalarField());
    }

    Info << "calculating cellLeafArea..." << endl;
    label celli;
    point faceCenter;
    vector faceArea;
    forAll(faces, i)
    {
        faceCenter = faces[i].centre(surfacePoints);
        celli = mesh.findCell(faceCenter);

        if (celli == -1){continue;}

        faceArea = surface_.faceAreas()[i];
        cellLeafArea[celli] += vector
                              (
                                sqrt(sqr(faceArea.x())),
                                sqrt(sqr(faceArea.y())),
                                sqrt(sqr(faceArea.z()))
                              );
        cellDiameterList[celli].append
        (
          sqrt(mag(faceArea) / (2 * constant::mathematical::pi))
        );
    } 
    Info << "Computing leafAreaDensity..." << endl;
    HashTable<dimensionedVector, label> leafAreaDensity(canopyCellsIndex.size());

    labelList cellList = cellLeafArea.toc();
    forAll(cellList, i)
    {
        this->la().set
        (
          cellList[i], 
          dimensionedScalar
          (
            dimArea/dimVolume, 
            mag(cellLeafArea[cellList[i]])
          )
        );
        this->lad().set
        (
          cellList[i], 
          dimensionedVector
          (
            dimArea/dimVolume, 
            cellLeafArea[cellList[i]] / mesh.cellVolumes()[cellList[i]]
          )
        );
        this->ldia().set
        (
          cellList[i], 
          dimensionedScalar
          (
            dimLength, 
            average(cellDiameterList[cellList[i]])
          )
        );
    }
}
template<class BaseCanopyModel>
dimensionedScalarCellSet canopyTriSurfaceModel<BaseCanopyModel>::calcLaCov
(
    const polyMesh& mesh,
    const triSurface& surface,
    const labelHashSet& canopyCellsIndex,
    vector direction
)
{
    // Loading surface
    triSurfaceSearch ts(surface);
    vector testDirection = (-1) * direction;

    // Find faces in each cell
    const faceList& faces = surface.faces(); 
    scalarField faceAreas = mag(surface.faceAreas())();
    const pointField& faceCs = surface.faceCentres();
    Info << "Number of pointIndexHits: " << faceCs.size() << endl;
    dimensionedScalarCellSet laCov(canopyCellsIndex.size());
    HashTable<List<face>,label> facesCellList(canopyCellsIndex.size());
    HashTable<pointField,label> faceCsCellList(canopyCellsIndex.size());
    forAllConstIter(labelHashSet, canopyCellsIndex, iter)
    {
        laCov.set
        (
          iter.key(),
          dimensionedScalar(dimless, 0)
        );
        facesCellList.set   
        (
          iter.key(),
          List<face>()
        );
        faceCsCellList.set
        (
          iter.key(),
          pointField()
        );

    }

    // Aggregate the surface face index to each cell index
    Info << "Aggregating surface face index to each cell index" << endl;
    label celli;
    // [TODO] Iterate to populate the HashTables
    forAll(faceCs, i)
    {
        celli = mesh.findCell(faceCs[i]);
        if (celli == -1){continue;}
        facesCellList[celli].append(faces[i]);
        faceCsCellList[celli].append(faceCs[i]);
    }

    // Find leaf coverage ratio
    Info << "Calculating leaf coverage ratios for direction " 
         << direction
         << endl;
    forAllConstIter(labelHashSet, canopyCellsIndex, iter)
    {
        // Projected area of cell at this direction
        // Projected area is calculated by computing the dot product of
        // each face area with the negative of incoming direction vector.
        // If the dot product is negative, then don't add as it will be
        // invisible
        
        cell cellIter = mesh.cells()[iter.key()];
        vector testVector = testDirection * 2 * mag(mesh.cellVolumes()[iter.key()]);

        scalar ACell = 0;
        forAll(cellIter, i)
        {
            face facei = mesh.faces()[cellIter[i]];
            label faceOwner = mesh.faceOwner()[cellIter[i]];
            scalar Afacei
            (
                // Face area vector
                // Note: face::area<pointField> is used instead of 
                // facei.area() as the latter gives (0 0 0) for some faces,
                // reasons unknown for now
                face::area<pointField>(facei.points(mesh.points())) 
                // aligned to the test direction
                & normalised(testDirection)
                // correct with owner
                * ( faceOwner == iter.key() ? 1 : (-1) )             
            );
            ACell += Afacei > 0 ? Afacei : 0;
        }

        // Compute projected leaf area of canopy
        faceList cellFaces = facesCellList[iter.key()];
        pointField startPts = (faceCsCellList[iter.key()] + normalised(testDirection) * 0.01)() ;
        pointField endPts = (startPts + testVector)();

        // Create a subset of the triSurface for each cell
        // This is important to avoid causing the error message flushing
        HashTable<label, label> cellFacesSurfacePtsTable(500);
        DynamicList<point> cellFacesSurfacePtsDyn(500);
        label ptCounter = 0;
        forAll(cellFaces, i)
        {
          forAll(cellFaces[i], j)
          {
              if(!cellFacesSurfacePtsTable.found(cellFaces[i][j]))
              {
                cellFacesSurfacePtsDyn.append(surface.points()[cellFaces[i][j]]);
                cellFacesSurfacePtsTable.insert(cellFaces[i][j], ptCounter);
                ptCounter++;
              }
          }
        }
        cellFacesSurfacePtsDyn.shrink();
        pointField cellFacesSurfacePts(cellFacesSurfacePtsDyn);
        
        List<labelledTri> cellFacesTri(cellFaces.size());
        forAll(cellFaces, i)
        {
            cellFacesTri[i] = 
            labelledTri
              (
                cellFacesSurfacePtsTable[cellFaces[i][0]],
                cellFacesSurfacePtsTable[cellFaces[i][1]],
                cellFacesSurfacePtsTable[cellFaces[i][2]],
                0
              )
            ;

        }
        triSurface cellTriSurface(cellFacesTri, cellFacesSurfacePts);

        triSurfaceSearch ts(cellTriSurface);
        List<pointIndexHit> hits(startPts.size());
        ts.findLine(startPts, endPts, hits);

        scalar ALeaves = 0;
        forAll(hits, i)
        {
          if (!hits[i].hit())
          {
            // Which side of the leaf is lit doesn't matter
            ALeaves += 
              mag
              (
                  normalised(testDirection) 
                & face::area<pointField>
                  (
                    cellFaces[i].points(surface.points())
                  )
              );
          }
        }
      
      laCov[iter.key()].value() = ALeaves / ACell;
    }

    return laCov;
    
}

template<class BaseCanopyModel>
dimensionedScalarCellSet canopyTriSurfaceModel<BaseCanopyModel>::calcLaLit
(
    const polyMesh& mesh,
    const triSurface& surface,
    const labelHashSet& canopyCellsIndex,
    vector direction
)
{
    // [TODO] Subset the surface bounded by the domain bbox to reduce work load
    triSurfaceSearch ts(surface);
    scalar zmax = max(mesh.points().component(2));
    vector testDirection = (-1) * direction;

    // Iterate over each leaf, 
    scalarField faceAreas = mag(surface.faceAreas())();
    const pointField& faceCs = surface.faceCentres();
    pointField startPts = (faceCs + normalised(testDirection) * 0.01)();
    const pointField endPts
    (
      startPts + testDirection * (zmax - startPts.component(2))
    );
    // Info << endPts << endl;

    Info << "Number of pointIndexHits: " << faceCs.size() << endl;
    dimensionedScalarCellSet laLit(canopyCellsIndex.size());
    forAllConstIter(labelHashSet, canopyCellsIndex, iter)
    {
        laLit.set
        (
          iter.key(),
          dimensionedScalar(dimArea, 0)
        );
    }
    List<pointIndexHit> hits(startPts.size());
    ts.findLine(startPts, endPts, hits);

    label cellj;
    forAll(hits, j)
    {
        if (hits[j].hit())
        {
        }
        else 
        {
            cellj = mesh.findCell(faceCs[j]);
            if (cellj == -1){continue;}
            laLit[cellj].value() += faceAreas[j];
        }
    }
    return laLit;
}




// template<class BaseCanopyModel>
// dimensionedVectorCellSet canopyTriSurfaceModel<BaseCanopyModel>::calcLAD
// (
//     const polyMesh& mesh,
//     const triSurface& surface,
//     const labelHashSet& cellsIndex
// )
// {
//     // Compute the total area of leaves within a mesh cell
//     const pointField& surfacePoints = surface.points();
//     const faceList& faces = surface.faces();
//     HashTable<vector, label> cellLeafArea(cellsIndex.size());
// 
//     forAll(cellsIndex, celli)
//     {
//         cellLeafArea.set(cellsIndex.toc()[celli], vector(0,0,0));
//     }
// 
//     Info << "calculating cellLeafArea..." << endl;
//     label celli;
//     point faceCenter;
//     vector faceArea;
//     forAll(faces, i)
//     {
//         faceCenter = faces[i].centre(surfacePoints);
//         celli = mesh.findCell(faceCenter);
//         faceArea = surface.faceAreas()[i];
//         cellLeafArea[celli] += vector
//                               (
//                                 sqrt(sqr(faceArea.x())),
//                                 sqrt(sqr(faceArea.y())),
//                                 sqrt(sqr(faceArea.z()))
//                               );
//     } 
//     Info << "Computing leafAreaDensity..." << endl;
//     HashTable<dimensionedVector, label> leafAreaDensity(cellsIndex.size());
// 
//     labelList cellList = cellLeafArea.toc();
//     forAll(cellList, i)
//     {
//         leafAreaDensity.set
//         (
//           cellList[i], 
//           dimensionedVector
//           (
//             dimArea/dimVolume, 
//             cellLeafArea[cellList[i]] / mesh.cellVolumes()[cellList[i]]
//           )
//         );
//     }
// 
//     return leafAreaDensity;
// }

template<class BaseCanopyModel>
canopyTriSurfaceModel<BaseCanopyModel>::canopyTriSurfaceModel
(
    const treeModel& tree
)
:
    canopySurfaceModel<BaseCanopyModel>(tree),
    surface_(fileName(this->surfaceModelDict().lookup("file")))
{
    init();
    // Testing find lalit
    vector defaultDirection(0, 0, -1);
    this->laLit() = calcLaLit
                    (
                      tree.mesh(), 
                      surface_, 
                      this->canopyCells(),
                      defaultDirection
                    );
    this->laCov() = calcLaCov
                    (
                      tree.mesh(), 
                      surface_, 
                      this->canopyCells(),
                      defaultDirection
                    );
};

} 

template<class BaseCanopyModel>
dimensionedScalarCellSet canopyTriSurfaceModel<BaseCanopyModel>::correctLaCov
(
  vector direction
)
{
    Info << "Correcting laCov for direction " << direction << endl;
    return calcLaCov
    (
      BaseCanopyModel::tree().mesh(),
      surface_,
      this->canopyCells(),
      direction
    );
}
