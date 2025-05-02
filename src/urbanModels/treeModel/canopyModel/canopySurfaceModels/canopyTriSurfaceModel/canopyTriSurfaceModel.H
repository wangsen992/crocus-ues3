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

Description
    Provides an interface for providing useful geoemtry information for
    the tree canopy. Offers both fine grained and bulk information based on
    different type of input. 

    LAD is calculated as scalar, not vector (isotropic lad)


SourceFiles
    canopyGeometryModel.C
\*---------------------------------------------------------------------------*/

#ifndef canopyTriSurfaceModel_H
#define canopyTriSurfaceModel_H

#include "fvCFD.H"
#include "treeModel.H"
#include "canopySurfaceModel.H"
#include "triSurfaceMesh.H"
#include "cellSetHashTables.H"

namespace Foam
{

template<class BaseCanopyModel>
class canopyTriSurfaceModel
  : public canopySurfaceModel<BaseCanopyModel>
{
    
    
    // Private Data

        triSurface surface_;

    // Private Member Functions

        void init();

        // Compute total leaf area lit in each cell with incoming 
        // light direction
        static dimensionedScalarCellSet calcLaLit
        (
            const polyMesh& mesh,
            const triSurface& surface,
            const labelHashSet& canopyCellsIndex,
            vector direction
        );
          
        static dimensionedScalarCellSet calcLaCov
        (
            const polyMesh& mesh,
            const triSurface& surface,
            const labelHashSet& canopyCellsIndex,
            vector direction
        );
            
        // static dimensionedVectorCellSet calcLAD
        // (
        //     const polyMesh& mesh,
        //     const triSurface& surface,
        //     const labelHashSet& cellsIndex
        // );

public: 
    
    // Constructors
    
        // Construct from surface file and mesh
        canopyTriSurfaceModel
        (
            const treeModel& tree
        );

    // Destructor
        ~canopyTriSurfaceModel(){};

    // Member Functions
        
        virtual dimensionedScalarCellSet correctLaCov
        (
          vector direction
        );

};
}


#ifdef NoRepository
  #include "canopyTriSurfaceModel.C"
#endif

#endif 
