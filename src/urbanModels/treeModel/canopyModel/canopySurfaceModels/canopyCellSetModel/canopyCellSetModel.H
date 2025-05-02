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

#ifndef canopyCellSetModel_H
#define canopyCellSetModel_H

#include "fvCFD.H"
#include "treeModel.H"
#include "canopySurfaceModel.H"
#include "cellSetHashTables.H"
#include "cellSet.H"

namespace Foam
{

template<class BaseCanopyModel>
class canopyCellSetModel
  : public canopySurfaceModel<BaseCanopyModel>
{
    
    
    // Private Data

        cellSet  set_;

    // Private Member Functions

        void init();


public: 
    
    // Constructors
    
        // Construct from surface file and mesh
        canopyCellSetModel
        (
            const treeModel& tree
        );

    // Destructor
        ~canopyCellSetModel(){};

    // Member Functions
        
        virtual dimensionedScalarCellSet correctLaCov
        (
          vector direction
        );

};
}


#ifdef NoRepository
  #include "canopyCellSetModel.C"
#endif

#endif 
