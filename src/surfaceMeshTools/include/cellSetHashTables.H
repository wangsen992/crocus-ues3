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

    To-Do: Inherit the HashTable to create new base class to allow additional
    operators, and methods such as .values() to get a full list of values

\*---------------------------------------------------------------------------*/

#ifndef cellSetHashTables_H
#define cellSetHashTables_H
#include "HashTable.H"
#include "dimensionedTypes.H"
#include "pointField.H"

namespace Foam
{
    typedef HashTable<scalar, label> scalarCellSet;
    typedef HashTable<vector, label> vectorCellSet;
    typedef HashTable<tensor, label> tensorCellSet;
    typedef HashTable<dimensionedScalar, label> dimensionedScalarCellSet;
    typedef HashTable<dimensionedVector, label> dimensionedVectorCellSet;
    typedef HashTable<dimensionedTensor, label> dimensionedTensorCellSet;

    // Convert the vectorCellSet into a list/field
    inline pointField vectorCellSet2pointField(vectorCellSet& vc)
    {
        const labelList& cellList(vc.sortedToc());
        pointField pF(cellList.size());
        label celli;
        forAll(cellList, i)
        {
            celli = cellList[i];
            pF[i] = vc[celli];
        }

        return pF;
    }
}

#endif
