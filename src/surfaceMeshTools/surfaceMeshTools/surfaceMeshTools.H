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
    Utility functions to help making surface-mesh type of operation easier.

SourceFiles
    surfaceMeshTools.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceMeshTools_H
#define surfaceMeshTools_H

#include "fvCFD.H"
#include "polyMesh.H"
#include "triSurfaceMesh.H"

namespace Foam
{
namespace surfaceMeshTools
{
    
    // Find mesh cells that are cut by surface
    labelHashSet findSurfaceCutCells
    (
        const polyMesh& mesh,
        const triSurfaceMesh& surface
    );

    labelHashSet findSurfaceCutCells
    (
        const polyMesh& mesh,
        const triSurface& surface
    );
    // Write the triSurface representation of the surface
    void writeObjFile(OFstream, triSurface&);
    
    // Write the polyMesh as obj file
    void writeObjFile(OFstream, polyMesh&);

    // Write obj representation of faces subset
    void writeObjFile
    (
        Ostream& os,
        const faceList& faces,
        const pointField& points,
        const labelList& faceLabels
    );

    // Write the obj representation of faces
    void writeObjFile
    (
        Ostream& os,
        const faceList&,
        const pointField&
    );

    // Return all faces intersected by or intersecting surface
    // Adapted from cellClassification::markFaces
    boolList markFaces
    (
      const polyMesh& mesh,
      const triSurface& surface,
      label debug=0
    );

    // Find intersection of an edge from two edge points and a face using
    // intersection from face
}
}

#endif
