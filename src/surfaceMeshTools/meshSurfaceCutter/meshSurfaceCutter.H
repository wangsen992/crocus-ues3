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

Class
    Foam::meshSurfaceCutter

Description
    Helper class that cuts a surface (comprises of faces) by mesh cells 
    & faces. 

SourceFiles
    meshSurfaceCutter.C
\*---------------------------------------------------------------------------*/

#ifndef meshSurfaceCutter_H
#define meshSurfaceCutter_H

#include "polyMesh.H"
#include "triSurfaceMesh.H"
#include "cellClassification.H"
#include "meshSearch.H"
#include "triSurfaceSearch.H"

namespace Foam
{

namespace surfaceMeshTools
{

class meshSurfaceCutter
{
    // Private Member Data

        //- Reference to mesh and surface
        const polyMesh&     mesh_;
        const triSurface&     surf_;

        //- cut information
        meshSearch queryMesh_;
        triSurfaceSearch  querySurf_;
        cellClassification cellClassification_;
        cellList meshCellsCutBySurf_;
        faceList meshFacesCutBySurf_;
        pointField intersectionPoints_;
        
    // Private Member Functions

        //- find intersections points between mesh and surf
        void calcIntersectionPoints
        (
            const polyMesh& mesh,
            const triSurface& surf_,
            pointField& intersectionPoints
        );


public:
    
    // Constructors

        //- Construt from components
        meshSurfaceCutter
        (
            const polyMesh& mesh, 
            const triSurface& surf
        );

    // Destructor
        
        ~meshSurfaceCutter(){};

    // Public Member Functions
        
        //- Find all intersection points on surface cut by mesh
        //- See tests/meshSurfaceCutter for current version
        pointField&  intersectionPoints(){return intersectionPoints_;};

        //- Cut surface with mesh
        tmp<triSurfaceMesh> cutSurf();

};
}
}

#endif
