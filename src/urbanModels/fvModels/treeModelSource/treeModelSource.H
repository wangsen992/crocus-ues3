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
    Foam::fv::treeModel

Description
    Interface to include tree models into the transport equations. 

    Currently the model is designed to be used for flexible runTime selection
    with different types of 

    The code structure (using templates and apply function) is borrowed from
    solidificationMeltingSource in the original source. 

    

Usage
    Example usage:
    \verbatim
    tree1
    {
        type        treeModel;

        U           U;  // Name of the velocity field
    }
    \endverbatim

SourceFiles
    buoyancyForce.C

\*---------------------------------------------------------------------------*/

#ifndef treeModelSource_H
#define treeModelSource_H

#include "fvModel.H"
#include "treeModel.H"
#include "treeModelList.H"

namespace Foam
{
namespace fv
{
/*---------------------------------------------------------------------------*\
                        Class treeModelSource Declaration
\*---------------------------------------------------------------------------*/

class treeModelSource
:
    public fvModel
{
    // Private Data

        const fvMesh& mesh_;
        IOdictionary urbanDict_;
        treeModelList treeList_;
        treeModel& tree_;
        bool updated_;
        label curTimeIndex_;

    // Private Member Functions
    
        dimensionedScalarCellSet getScalarSource(const word&) const;

public:

    //- Runtime type information
    TypeName("treeModelSource");


    // Constructors

        //- Construct from explicit source name and mesh
        treeModelSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Disallow default bitwise copy construction
        treeModelSource(const treeModelSource&) = delete;
        
    // Member Functions

        // Checks

            //- Return the list of fields for which the fvModel adds source term
            //  to the transport equation
            virtual wordList addSupFields() const;


        // Evaluate

             //- Add explicit contribution to incompressible momentum equation
             virtual void addSup
             (
                 fvMatrix<vector>& eqn,
                 const word& fieldName
             ) const;
 
             //- Add explicit contribution to compressible momentum equation
             virtual void addSup
             (
                 const volScalarField& rho,
                 fvMatrix<vector>& eqn,
                 const word& fieldName
             ) const;
 
             //- Add explicit contribution to phase momentum equation
             virtual void addSup
             (
                 const volScalarField& alpha,
                 const volScalarField& rho,
                 fvMatrix<vector>& eqn,
                 const word& fieldName
             ) const;
 
             //- Add explicit contribution to incompressible scalar equation
             virtual void addSup
             (
                 fvMatrix<scalar>& eqn,
                 const word& fieldName
             ) const;
 
             //- Add explicit contribution to compressible scalar equation
             virtual void addSup
             (
                 const volScalarField& rho,
                 fvMatrix<scalar>& eqn,
                 const word& fieldName
             ) const;
 
             //- Add explicit contribution to phase scalar equation
             virtual void addSup
             (
                 const volScalarField& alpha,
                 const volScalarField& rho,
                 fvMatrix<scalar>& eqn,
                 const word& fieldName
             ) const;

             virtual void correct();

        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const treeModelSource&) = delete;
};


}
}

#endif
