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
    Foam::canopySurfaceModel

Description
    Base class to provide access to LAD related info


SourceFiles
    canopySurfaceModel.C
\*---------------------------------------------------------------------------*/

#ifndef canopySurfaceModel_H
#define canopySurfaceModel_H

#include "fvMesh.H"
#include "treeModel.H"
#include "cellSetHashTables.H"

namespace Foam
{

template<class BaseCanopyModel>
class canopySurfaceModel
:
    public BaseCanopyModel
{
    
    // Private Member Data
        //- Surface Model dictionary
        
        //- cell list influenced by the canopy
        labelHashSet     canopyCellsIndex_;

        //- Leaf area density
        volScalarField lad_;
        
        //- Total Leaf Area
        volScalarField la_;

        //- Total Sunlit Leaf Area
        volScalarField laLit_;

        //- Leaf coverage ratio for each cell (between 0 and 1);
        volScalarField laCov_;

        //- Cell-Averaged Leave Diameter
        // Note: Currently assume all leaves are spherical!!
        volScalarField ldia_;
        
    // Private Member Functions


public:
        
    // Constructors
    
        // Construct from surface file and mesh
        canopySurfaceModel
        (
            const treeModel& tree
        );

    // Destructor
        ~canopySurfaceModel(){};

    // Member Functions

        // Access
            //- reference to the surface model dict

        // Parameters for all types of analysis

            //- reference to canopy cells
            virtual const labelHashSet& canopyCells() const {return canopyCellsIndex_;};
            virtual labelHashSet& canopyCells()  {return canopyCellsIndex_;};
            //- reference to leaf area density
            virtual volScalarField& lad() {return lad_;};
            virtual const volScalarField& lad() const {return lad_;};

            virtual volScalarField& la() {return la_;};
            virtual const volScalarField& la() const {return la_;};

            virtual volScalarField& laLit() {return laLit_;};
            virtual const volScalarField& laLit() const {return laLit_;};
            virtual volScalarField& laCov() {return laCov_;};
            virtual const volScalarField& laCov() const {return laCov_;};
            
            virtual volScalarField& ldia() {return ldia_;};
            virtual const volScalarField& ldia() const {return ldia_;};

        // Correction functions
            virtual volScalarField correctLaCov
            (
              vector direction
            ){return laCov_;};
};
}

#ifdef NoRepository
  #include "canopySurfaceModel.C"
#endif


#endif
