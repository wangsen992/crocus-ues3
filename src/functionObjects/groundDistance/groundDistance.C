/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2024 OpenCFD Ltd.
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

#include "groundDistance.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(groundDistance, 0);
    addToRunTimeSelectionTable(functionObject, groundDistance, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::groundDistance::groundDistance
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict)
{
    read(dict);

    volScalarField* distPtr
    (
        new volScalarField
        (
            IOobject
            (
                "groundDistance",
                mesh_.time().constant(),
                mesh_.thisDb(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::REGISTER
            ),
            mesh_,
            dimensionedScalar(dimLength, Zero)
        )
    );

    regIOobject::store(distPtr);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::groundDistance::read
(
    const dictionary& dict
)
{
    fvMeshFunctionObject::read(dict);

    doCells_ = dict.getOrDefault("calculateCells", true);

    geomPtr_.reset(nullptr);
    geomPtr_.reset
    (
        new searchableSurfaces
        (
            IOobject
            (
                "abc",                             // dummy name
                mesh_.time().constant(),           // directory
                "triSurface",                      // instance
                mesh_.time(),                      // registry
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            dict.subDict("geometry"),
            true                // allow single-region shortcut
        )
    );

    return true;
}


bool Foam::functionObjects::groundDistance::execute()
{
    volScalarField& distance = mesh_.lookupObjectRef<volScalarField>
    (
        "groundDistance"
    );

    volScalarField::Boundary& bfld = distance.boundaryFieldRef();
    forAll(bfld, patchi)
    {
        if (!polyPatch::constraintType(bfld[patchi].patch().type()))
        {
            const pointField& fc = mesh_.C().boundaryField()[patchi];

            labelList surfaces;
            List<pointIndexHit> nearestInfo;
            geomPtr_().findAnyIntersection
            (
                fc,
                fc + vector(0,0,-1000),
                surfaces,
                nearestInfo
            );
            // geomPtr_().findNearest
            // (
            //     fc,
            //     scalarField(fc.size(), GREAT),
            //     surfaces,
            //     nearestInfo
            // );

            scalarField dist(fc.size());
            forAll(nearestInfo, i)
            {
                if (nearestInfo[i].hit())
                {
                  dist[i] = nearestInfo[i].hitPoint().dist(fc[i]);
                }
                else
                {
                  dist[i] = 0;
                }
            }
            bfld[patchi] == dist;
        }
    }

    if (doCells_)
    {
        const pointField& cc = mesh_.C().internalField();

        labelList surfaces;
        List<pointIndexHit> nearestInfo;
        geomPtr_().findNearest
        (
            cc,
            scalarField(cc.size(), GREAT),
            surfaces,
            nearestInfo
        );

        forAll(nearestInfo, celli)
        {
            distance[celli] = nearestInfo[celli].hitPoint().dist(cc[celli]);
        }
    }
    distance.correctBoundaryConditions();

    return true;
}


bool Foam::functionObjects::groundDistance::write()
{
    Log << "    functionObjects::" << type() << " " << name()
        << " writing distance-to-surface field" << endl;

    const volScalarField& distance =
        mesh_.lookupObject<volScalarField>("groundDistance");

//    volScalarField::Boundary& bfld = distance.boundaryFieldRef();
//    forAll(bfld, patchi)
//    {
//        if (!polyPatch::constraintType(bfld[patchi].patch().type()))
//        {
//            const pointField& fc = mesh_.C().boundaryField()[patchi];
//
//            labelList surfaces;
//            List<pointIndexHit> nearestInfo;
//            geomPtr_().findNearest
//            (
//                fc,
//                scalarField(fc.size(), GREAT),
//                surfaces,
//                nearestInfo
//            );
//
//            scalarField dist(fc.size());
//            forAll(nearestInfo, i)
//            {
//                dist[i] = nearestInfo[i].hitPoint().dist(fc[i]);
//            }
//            bfld[patchi] == dist;
//        }
//    }
//
//    if (doCells_)
//    {
//        const pointField& cc = mesh_.C().internalField();
//
//        labelList surfaces;
//        List<pointIndexHit> nearestInfo;
//        geomPtr_().findNearest
//        (
//            cc,
//            scalarField(cc.size(), GREAT),
//            surfaces,
//            nearestInfo
//        );
//
//        forAll(nearestInfo, celli)
//        {
//            distance[celli] = nearestInfo[celli].hitPoint().dist(cc[celli]);
//        }
//    }
//    distance.correctBoundaryConditions();
    distance.write();

    return true;
}


// ************************************************************************* //
