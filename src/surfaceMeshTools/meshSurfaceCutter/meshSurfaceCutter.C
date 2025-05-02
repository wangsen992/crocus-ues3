/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2021 OpenFOAM Foundation
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

Description
    Helper class that cuts a surface (comprises of faces) by mesh cells 
    & faces. 

\*---------------------------------------------------------------------------*/
    
#include "meshSurfaceCutter.H"
#include "surfaceMeshTools.H"
#include "meshSearch.H"
#include "DynamicList.H"
#include "indexedOctree.H"
#include "treeDataFace.H"
#include "treeBoundBox.H"

namespace Foam
{
namespace surfaceMeshTools
{

void meshSurfaceCutter::calcIntersectionPoints
(
    const polyMesh& mesh,
    const triSurface& surf,
    pointField& intersectionPoints
)
{
    const List<labelledTri>& surfFaces(surf.surfFaces());
    labelList meshFacesList(mesh.nInternalFaces());
    DynamicList<point> intersectionPointsDynList(0);
    forAll(meshFacesList, i)
    {
        meshFacesList[i] = i;
    }
    indexedOctree<treeDataFace> faceTree
    (
      treeDataFace
      (false, mesh, meshFacesList),
      treeBoundBox(mesh.points()),
      8, 10, 3.00
    );
    forAll(surfFaces, i)
    {
        face facei = surfFaces[i];
        const pointField& facePoints = facei.points(surf.points());
        forAll(facePoints, j)
        {
            point p1 = facePoints[j];
            point p2 = facePoints[facePoints.fcIndex(j)];
            pointIndexHit hit = faceTree.findLine(p1, p2);
            if (hit.hit())
            {
                intersectionPointsDynList.append(hit.hitPoint());
            }
        }
    }

    intersectionPoints.resize(intersectionPointsDynList.size());
    forAll(intersectionPoints, i)
    {
        intersectionPoints[i] = intersectionPointsDynList[i];
    }
}

meshSurfaceCutter::meshSurfaceCutter
(
    const polyMesh& mesh, 
    const triSurface& surf
)
:
mesh_(mesh),
surf_(surf),
queryMesh_(mesh_),
querySurf_(surf_),
cellClassification_
(
  mesh_,queryMesh_,querySurf_,pointField(1, vector(0, 0, 0))
),
meshCellsCutBySurf_(),
meshFacesCutBySurf_(),
intersectionPoints_(0)
{
    forAll(cellClassification_, i)
    {
        if (cellClassification_[i] == cellClassification::CUT)
        {
            meshCellsCutBySurf_.append(mesh.cells()[i]);
            
        };
    }

    // find all intersection points
    calcIntersectionPoints(mesh_, surf_, intersectionPoints_);
};

tmp<triSurfaceMesh> meshSurfaceCutter::cutSurf()
{
    // Set up containers for following operation
    labelHashSet cutCellSet = findSurfaceCutCells(mesh_, surf_);
    labelList cutCells = cutCellSet.sortedToc();
    
    // For each triangle, there are only several cases as listed below
    // Case 1: All three points are in one cells
    // Case 2: One cell has 1 point, the neighbor cell has 2 points
    // Case 3: three cells each contains 1 point
    // For each case, there is also the possibility that the surface face
    // spread far such that there are cells cut by surface but contain no
    // vertex

    


    // Split up surface faces (trianglized) that are cut by mesh faces
    // Check each surface face
    // If contained, leave alone
    // If cut through n face(s), split face into n+1 sub-faces

    // The final result is a HashTable with cell index as key and a list
    // of faces within that cell as the value
}
}
}
