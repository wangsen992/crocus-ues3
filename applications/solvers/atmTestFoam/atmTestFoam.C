/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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
    buoyantBoussinesqPimpleFoam

Group
    grpHeatTransferSolvers

Description
    Transient solver for buoyant, turbulent flow of incompressible fluids,
    with optional mesh motion and mesh topology changes.

    Uses the Boussinesq approximation:
    \f[
        rho_{k} = 1 - beta(T - T_{ref})
    \f]

    where:
        \f$ rho_{k} \f$ = the effective (driving) kinematic density
        beta = thermal expansion coefficient [1/K]
        T = temperature [K]
        \f$ T_{ref} \f$ = reference temperature [K]

    Valid when:
    \f[
        \frac{beta(T - T_{ref})}{rho_{ref}} << 1
    \f]

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "radiationModel.H"
#include "CorrectPhi.H"
#include "fvOptions.H"
#include "pimpleControl.H"
#include "interpolate2D.H"
#include "windRoseToCartesian.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for buoyant, turbulent flow"
        " of incompressible fluids, with optional mesh"
        " motion and mesh topology changes.\n"
        "Uses the Boussinesq approximation."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "createDyMControls.H"
    #include "createFields.H"

    Info << "Creating post processing dir" << endl;
    #include "createPostProcessingDir.H"
    Info << "find vertical cell levels" << endl;
    #include "findVerticalCellLevels.H"
    Info << "read abl properties" << endl;
    #include "readABLProperties.H"
    Info << "create source terms" << endl;
    #include "createSourceTerms.H"

    #include "createUfIfPresent.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"
    #include "initContinuityErrs.H"

    turbulence->validate();
    
    // create coefficients for lateral boundary damping
    volScalarField dampingCoeff
    (
        IOobject
        (
            "dampingCoeff",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("dampingCoeff",dimless/dimTime,0.0)
    );

    scalar Lx(ABLProperties.lookupOrDefault<scalar>("xMax", 1000.0));
    scalar Ly(ABLProperties.lookupOrDefault<scalar>("yMax", 1000.0));
    scalar d(ABLProperties.lookupOrDefault<scalar>("dDamp", 200.0));
    scalar maxC(ABLProperties.lookupOrDefault<scalar>("CDamp", 1e-5));
    forAll(dampingCoeff, celli)
    {
        scalar x = mesh.cellCentres()[celli].x();
        scalar y = mesh.cellCentres()[celli].y();
        scalar z = mesh.cellCentres()[celli].z();

        scalar cur = dampingCoeff.primitiveFieldRef()[celli];
        if (x > (Lx - d))
        {
          dampingCoeff.primitiveFieldRef()[celli] = max((d - (Lx - x))/d * maxC, cur);
        }
        else if (x < (d - Lx))
        {
          dampingCoeff.primitiveFieldRef()[celli] = max((d - (Lx + x))/d * maxC, cur);
        }
        else if (y > (Ly - d))
        {
          dampingCoeff.primitiveFieldRef()[celli] = max((d - (Ly - y))/d * maxC, cur);
        }
        else if (y < (d - Ly))
        {
          dampingCoeff.primitiveFieldRef()[celli] = max((d - (Ly + y))/d * maxC, cur);
        }
    }

      

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                // Do any mesh changes
                mesh.controlledUpdate();

                if (mesh.changing())
                {
                    MRF.update();

                    if (correctPhi)
                    {
                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & Uf();

                        #include "correctPhi.H"

                        // Make the flux relative to the mesh motion
                        fvc::makeRelative(phi, U);
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }

            #include "UEqn.H"
            #include "TEqn.H"
            #include "qEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                laminarTransport.correct();
                turbulence->correct();
            }

            #include "correctSourceTerms.H"
        }

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
