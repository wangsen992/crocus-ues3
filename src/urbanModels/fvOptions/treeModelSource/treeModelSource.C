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

\*---------------------------------------------------------------------------*/

#include "treeModelSource.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(treeModelSource, 0);

    addToRunTimeSelectionTable
    (
        option,
        treeModelSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
tmp<volScalarField> Foam::fv::treeModelSource::getScalarSource(const word& fieldName) const
{
    // dimensionedScalarCellSet Fi(tree_.canopy().canopyCells().size());
    tmp<volScalarField> Fi;
    Info << "Selecting treeSource fieldName : " << fieldName << endl;
    if (fieldName == "k")
    {
        Fi.clear();
        Fi = tree_.canopy().Fturb("k");
    }
    else if (fieldName == "epsilon")
    {
        Fi.clear();
        Fi = tree_.canopy().Fturb("epsilon");
    }
    else if (fieldName == "T")
    {
        Fi.clear();
        Fi = tree_.canopy().FT();
    }
    else if (fieldName == "H2O")
    {
        Fi.clear();
        Fi = tree_.canopy().Fq();
    }

    return Fi;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::treeModelSource::treeModelSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    option(name, modelType, dict, mesh),
    mesh_(mesh),
    urbanDict_
    (
      IOobject
      (
        "urbanModelProperties",
        mesh.time().constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
      )
    ),
    treeList_
    (
        mesh, 
        urbanDict_.subDict("treeList")
    ),
    tree_
    (
      treeList_[0]
    ),
    updated_(false),
    curTimeIndex_(mesh.time().timeIndex())
{
    Info << "[debug] treeModelSource Init: " << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::treeModelSource::addSupFields() const
{
    wordList supFields = {"U", "T", "H2O"};
    return supFields;
}

void Foam::fv::treeModelSource::addSup
(
    fvMatrix<vector>& eqn,
    const word& fieldName
) const
{
    Info << "treeModelSource addSup: " << fieldName << endl;
    vectorField& Usource =  eqn.source();
    auto Fu = tree_.canopy().Fu();
    Usource -= Fu * mesh_.V();
}


void Foam::fv::treeModelSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const word& fieldName
) const
{
    Info << "treeModelSource addRhoSup: " << fieldName << endl;
    vectorField& Usource =  eqn.source();
    auto Fu = tree_.canopy().Fu();
    Usource -= rho * mesh_.V() * Fu;
}


void Foam::fv::treeModelSource::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const word& fieldName
) const
{
    Info << "treeModelSource addAlphaRhoSup: " << fieldName << endl;
    vectorField& Usource =  eqn.source();
    auto Fu = tree_.canopy().Fu();
    Usource -= alpha * rho * mesh_.V() * Fu;
}

// Correcting all scalar equations
void Foam::fv::treeModelSource::addSup
(
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{

    Info << "treeModelSource addSupS: "<< fieldName  << endl;
    scalarField& source =  eqn.source();
    auto Fi = getScalarSource(fieldName);
    source -= Fi * mesh_.V();
        
}


void Foam::fv::treeModelSource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    Info << "treeModelSource addRhoSupC: " << fieldName << endl;
    scalarField& source =  eqn.source();
    auto Fi = getScalarSource(fieldName);
    source -= rho  * mesh_.V() * Fi;
}


void Foam::fv::treeModelSource::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    Info << "treeModelSource addAlphaRhoSupC: "<< fieldName  << endl;
    scalarField& source =  eqn.source();
    auto Fi = getScalarSource(fieldName);
    source -= alpha * rho * mesh_.V() * Fi;
}

void Foam::fv::treeModelSource::correct()
{
    if(curTimeIndex_ != this->mesh().time().timeIndex())
    {
      curTimeIndex_ = this->mesh().time().timeIndex();
      Info << "TreeModelSource Correcting.." << endl;
      tree_.canopy().correctMomentumTransfer();
      tree_.canopy().correctEnergyTransfer();
    }
}

bool Foam::fv::treeModelSource::read(const dictionary& dict)
{
    if (option::read(dict))
    {
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
