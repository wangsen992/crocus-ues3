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

\*---------------------------------------------------------------------------*/

#include "labelledTri.H"
#include "surfaceMeshTools.H"
#include "triSurfaceSearch.H"
#include "treeDataFace.H"
#include "meshTools.H"

namespace Foam
{

labelHashSet surfaceMeshTools::findSurfaceCutCells
(
    const polyMesh &mesh, 
    const triSurfaceMesh &surfaceMesh
)
{
    triSurface surface(surfaceMesh.surface());

    return findSurfaceCutCells(mesh, surface);
}

labelHashSet surfaceMeshTools::findSurfaceCutCells
(
    const polyMesh &mesh, 
    const triSurface &surface
)
{
    Info << "running findSurfaceCutCells" << endl; 
    DynamicList<label> dynList(0);

    // Compute the total area of leaves within a mesh cell
    const pointField& surfacePoints = surface.points();
    const List<labelledTri>& faces = surface.surfFaces();

    label celli;
    point faceCenter;
    forAll(faces, i)
    {
        
        faceCenter = faces[i].centre(surfacePoints);
        celli = mesh.findCell(faceCenter);
        // Defend against celli == -1
        if (celli == -1){continue;}
        dynList.append(celli);
    } 

    return labelHashSet(dynList);
}
void surfaceMeshTools::writeObjFile(OFstream of, triSurface& surf)
{
    forAll(surf.points(), i)
    {
        of << "v " 
              << surf.points()[i].x() << " "
              << surf.points()[i].y() << " "
              << surf.points()[i].z() << endl;
    }
    forAll(surf.surfFaces(), i)
    {
        face facei = surf.surfFaces()[i];
        of << "f ";
        forAll(facei, plabel)
        {
          of << facei[plabel]+1 << " ";
        }
        of << endl;
    }
    // To-DO: Check if the close file is needed.
}

void surfaceMeshTools::writeObjFile(OFstream of, polyMesh& mesh)
{
    forAll(mesh.points(), i)
    {
        of << "v " 
              << mesh.points()[i].x() << " "
              << mesh.points()[i].y() << " "
              << mesh.points()[i].z() << endl;
    }
    forAll(mesh.faces(), i)
    {
        face facei = mesh.faces()[i];
        of << "f ";
        forAll(facei, plabel)
        {
          of << facei[plabel]+1 << " ";
        }
        of << endl;
    }
    // To-DO: Check if the close file is needed.
}

boolList surfaceMeshTools::markFaces
(
  const polyMesh& mesh, 
  const triSurface& surface,
  label debug
)
{
    triSurfaceSearch search(surface);
    cpuTime timer;

    boolList cutFace(mesh.nFaces(), false);

    label nCutFaces = 0;

    // Intersect mesh edges with surface (is fast) and mark all faces that
    // use edge.

    forAll(mesh.edges(), edgeI)
    {
        if (debug && (edgeI % 10000 == 0))
        {
            Pout<< "Intersecting mesh edge " << edgeI << " with surface"
                << endl;
        }

        const edge& e = mesh.edges()[edgeI];

        const point& p0 = mesh.points()[e.start()];
        const point& p1 = mesh.points()[e.end()];

        pointIndexHit pHit(search.tree().findLineAny(p0, p1));

        if (pHit.hit())
        {
            const labelList& myFaces = mesh.edgeFaces()[edgeI];

            forAll(myFaces, myFacei)
            {
                label facei = myFaces[myFacei];

                if (!cutFace[facei])
                {
                    cutFace[facei] = true;

                    nCutFaces++;
                }
            }
        }
    }

    if (debug)
    {
        Pout<< "Intersected edges of mesh with surface in = "
            << timer.cpuTimeIncrement() << " s\n" << endl << endl;
    }

    //
    // Construct octree on faces that have not yet been marked as cut
    //

    labelList allFaces(mesh.nFaces() - nCutFaces);

    label allFacei = 0;

    forAll(cutFace, facei)
    {
        if (!cutFace[facei])
        {
            allFaces[allFacei++] = facei;
        }
    }

    if (debug)
    {
        Pout<< "Testing " << allFacei << " faces for piercing by surface"
            << endl;
    }

    treeBoundBox allBb(mesh.points());

    // Extend domain slightly (also makes it 3D if was 2D)
    scalar tol = 1e-6 * allBb.avgDim();

    point& bbMin = allBb.min();
    bbMin.x() -= tol;
    bbMin.y() -= tol;
    bbMin.z() -= tol;

    point& bbMax = allBb.max();
    bbMax.x() += 2*tol;
    bbMax.y() += 2*tol;
    bbMax.z() += 2*tol;

    indexedOctree<treeDataFace> faceTree
    (
        treeDataFace(false, mesh, allFaces),
        allBb,      // overall search domain
        8,          // maxLevel
        10,         // leafsize
        3.0         // duplicity
    );

    const triSurface& surf = search.surface();
    const edgeList& edges = surf.edges();
    const pointField& localPoints = surf.localPoints();

    label nAddFaces = 0;

    forAll(edges, edgeI)
    {
        if (debug && (edgeI % 10000 == 0))
        {
            Pout<< "Intersecting surface edge " << edgeI
                << " with mesh faces" << endl;
        }
        const edge& e = edges[edgeI];

        const point& start = localPoints[e.start()];
        const point& end = localPoints[e.end()];

        vector edgeNormal(end - start);
        const scalar edgeMag = mag(edgeNormal);
        const vector smallVec = 1e-9*edgeNormal;

        edgeNormal /= edgeMag+SMALL;

        // Current start of pierce test
        point pt = start;

        while (true)
        {
            pointIndexHit pHit(faceTree.findLine(pt, end));

            if (!pHit.hit())
            {
                break;
            }
            else
            {
                label facei = faceTree.shapes().faceLabels()[pHit.index()];

                if (!cutFace[facei])
                {
                    cutFace[facei] = true;

                    nAddFaces++;
                }

                // Restart from previous endpoint
                pt = pHit.hitPoint() + smallVec;

                if (((pt-start) & edgeNormal) >= edgeMag)
                {
                    break;
                }
            }
        }
    }

    if (debug)
    {
        Pout<< "Detected an additional " << nAddFaces << " faces cut"
            << endl;

        Pout<< "Intersected edges of surface with mesh faces in = "
            << timer.cpuTimeIncrement() << " s\n" << endl << endl;
    }

    return cutFace;
}


void Foam::surfaceMeshTools::writeObjFile
(
    Ostream& os,
    const faceList& faces,
    const pointField& points,
    const labelList& faceLabels
)
{
    Map<label> foamToObj(4*faceLabels.size());

    label vertI = 0;

    forAll(faceLabels, i)
    {
        const face& f = faces[faceLabels[i]];

        forAll(f, fp)
        {
            if (foamToObj.insert(f[fp], vertI))
            {
                // Note that library native func is called to write vertices
                meshTools::writeOBJ(os, points[f[fp]]);
                vertI++;
            }
        }

        // Modification is made here from "l" to "f"
        os << 'f';
        forAll(f, fp)
        {
            os << ' ' << foamToObj[f[fp]]+1;
        }
        os <<  endl;
    }
}


void Foam::surfaceMeshTools::writeObjFile
(
    Ostream& os,
    const faceList& faces,
    const pointField& points
)
{
    labelList allFaces(faces.size());
    forAll(allFaces, i)
    {
        allFaces[i] = i;
    }
    writeObjFile(os, faces, points, allFaces);
}


}
