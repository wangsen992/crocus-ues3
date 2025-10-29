/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

Application
    csv2var

Group
    grpMeshManipulationUtilities

Description
    Divides external faces into patches based on (user supplied) feature
    angle.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "dimensionSet.H"
#include "polyMesh.H"
#include "Time.H"
#include "volFieldsFwd.H"
#include "fvCFD.H"
#include "meshSearch.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
int main(int argc, char *argv[])
{

    // argList::addArgument("var", "variable name (must be available in constant");

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    Info<< "Mesh read in = "
        << runTime.cpuTimeIncrement()
        << " s\n" << endl << endl;

    word var_name("lad");
    word data_path("./constant/urban/point_data");
    // Read var file
    IOdictionary var_dict
    (
      IOobject
      (
          var_name,
          data_path,
          runTime,
          IOobject::MUST_READ,
          IOobject::NO_WRITE
      )
    );
    scalarField var_data(var_dict.get<scalarField>(var_name));
    vectorField points(var_dict.get<vectorField>("points"));
    dimensionSet dim(var_dict.get<dimensionSet>("dimension"));
    for (label i=0;i<10;i++)
    {
      Info << points[i] <<", "<< var_data[i] << endl;
    }
    Info << dim << endl;

    // create var as geometricField
    volScalarField lad
    (
      IOobject
      (
        "canopy_lad",
        runTime.constant(),
        runTime,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
      ),
      mesh,
      dimensionedScalar(dim, 0.0)
    );
    volScalarField la
    (
      IOobject
      (
        "canopy_la",
        runTime.constant(),
        runTime,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
      ),
      mesh,
      dimensionedScalar(dim, 0.0)
    );
    volScalarField laLit
    (
      IOobject
      (
        "canopy_laLit",
        runTime.constant(),
        runTime,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
      ),
      mesh,
      dimensionedScalar(dim, 0.0)
    );
    volScalarField laCov
    (
      IOobject
      (
        "canopy_laCov",
        runTime.constant(),
        runTime,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
      ),
      mesh,
      dimensionedScalar(dim, 0.0)
    );
    volScalarField ldia
    (
      IOobject
      (
        "canopy_ldia",
        runTime.constant(),
        runTime,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
      ),
      mesh,
      dimensionedScalar(dim, 0.1)
    );

    meshSearch ms(mesh);

    scalar k = 0.5;
    scalar t = 0; // transmissivity
    label num_points = points.size();
    scalar perc = 0.05;
    forAll(points, i)
    {
      if (i > num_points * perc)
      {
          Pout << "processed " << perc * 100 << " points..." << endl;
          perc += 0.05;
      }
      if (!ms.isInside(points[i]))
      {
        continue;
      }
      label celli = ms.findCell(points[i]);
      scalar ladi;
      if (celli != -1)
      {
          // Info << celli << endl;
          ladi = min(1 * var_data[i], 3);
          lad.primitiveFieldRef()[celli] = ladi;
          la.primitiveFieldRef()[celli] = ladi * mesh.V()[celli];
          t = 1/Foam::exp(ladi * k * pow(mesh.V()[celli],1/3));
          laCov.primitiveFieldRef()[celli] = 1 - t ;
          laLit.primitiveFieldRef()[celli] = pow(mesh.V()[celli], 2/3) * (1-t);
      }
    }

    // More precision (for points data)
    IOstream::minPrecision(10);

    lad.setUpToDate();
    lad.write();
    la.setUpToDate();
    la.write();
    laLit.setUpToDate();
    laLit.write();
    laCov.setUpToDate();
    laCov.write();
    ldia.setUpToDate();
    ldia.write();
    // mesh.write();
    // runTime.setUpToDate();
    // runTime.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
