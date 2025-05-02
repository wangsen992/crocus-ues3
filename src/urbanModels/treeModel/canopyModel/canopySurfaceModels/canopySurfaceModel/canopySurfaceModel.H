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
        dictionary      surfaceModelDict_;
        
        //- cell list influenced by the canopy
        labelHashSet     canopyCellsIndex_;

        //- Leaf area density
        dimensionedVectorCellSet lad_;
        
        //- Total Leaf Area
        dimensionedScalarCellSet la_;

        //- Total Sunlit Leaf Area
        dimensionedScalarCellSet laLit_;

        //- Leaf coverage ratio for each cell (between 0 and 1);
        dimensionedScalarCellSet laCov_;

        //- Cell-Averaged Leave Diameter
        // Note: Currently assume all leaves are spherical!!
        dimensionedScalarCellSet ldia_;
        
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
            virtual const dictionary& surfaceModelDict() const {return surfaceModelDict_;}

        // Parameters for all types of analysis

            //- reference to canopy cells
            virtual const labelHashSet& canopyCells() const {return canopyCellsIndex_;};
            virtual labelHashSet& canopyCells()  {return canopyCellsIndex_;};
            //- reference to leaf area density
            virtual dimensionedVectorCellSet& lad() {return lad_;};
            virtual const dimensionedVectorCellSet& lad() const {return lad_;};

            virtual dimensionedScalarCellSet& la() {return la_;};
            virtual const dimensionedScalarCellSet& la() const {return la_;};

            virtual dimensionedScalarCellSet& laLit() {return laLit_;};
            virtual const dimensionedScalarCellSet& laLit() const {return laLit_;};
            virtual dimensionedScalarCellSet& laCov() {return laCov_;};
            virtual const dimensionedScalarCellSet& laCov() const {return laCov_;};
            
            virtual dimensionedScalarCellSet& ldia() {return ldia_;};
            virtual const dimensionedScalarCellSet& ldia() const {return ldia_;};

        // Correction functions
            virtual dimensionedScalarCellSet correctLaCov
            (
              vector direction
            ) = 0;
};
}

#ifdef NoRepository
  #include "canopySurfaceModel.C"
#endif


#endif
