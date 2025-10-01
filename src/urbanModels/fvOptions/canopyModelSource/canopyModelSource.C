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

#include "canopyModelSource.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(canopyModelSource, 0);

    addToRunTimeSelectionTable
    (
        option,
        canopyModelSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
tmp<volScalarField> Foam::fv::canopyModelSource::getScalarSource(const word& fieldName) const
{
    // dimensionedScalarCellSet Fi(canopy_.canopyCells().size());
    tmp<volScalarField> Fi;
    Info << "Selecting treeSource fieldName : " << fieldName << endl;
    if (fieldName == "k")
    {
        Fi.clear();
        Fi = canopy_->Fturb("k");
    }
    else if (fieldName == "epsilon")
    {
        Fi.clear();
        Fi = canopy_->Fturb("epsilon");
    }
    else if (fieldName == "T")
    {
        Fi.clear();
        Fi = canopy_->FT();
        Info << "[debug] check FT: " << max(Fi) << endl;
        Info << "[debug] check FT: " << min(Fi) << endl;
    }
    else if (fieldName == "H2O")
    {
        Fi.clear();
        Fi = canopy_->Fq();
        Info << "[debug] check Fq: " << average(Fi) << endl;
    }

    return Fi;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::canopyModelSource::canopyModelSource
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
    tree_
    (
        mesh, 
        urbanDict_.subDict("treeDict")
    ),
    canopy_
    (
      canopyModel::New(tree_)
    ),
    updated_(false),
    curTimeIndex_(mesh.time().timeIndex())
{
    Info << "[debug] canopyModelSource Init: " << endl;
    fieldNames_.resize(1, "U");
    fieldNames_.append("T");
    fieldNames_.append("H2O");
    Info << "[debug] canopyModelSource Init Complete " << endl;
    Info << "[debug] Fu: " << canopy_->addSupFields() << endl;
    fv::option::resetApplied();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Foam::wordList Foam::fv::canopyModelSource::addSupFields() const
// {
//     wordList supFields = {"U", "T", "H2O"};
//     return supFields;
// }

void Foam::fv::canopyModelSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
) 
{
    eqn += canopy_->Fu();
    // Info << "canopyModelSource addSup: " << fieldi << endl;
    // vectorField& Usource =  eqn.source();
    // Info << "canopyModelSource addSup eqnSource "<< endl;
    // volVectorField& Fu = canopy_->Fu();
    // Info << "canopyModelSource addSup Fu "<< endl;
    // Usource -= Fu * mesh_.V();
}


void Foam::fv::canopyModelSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    eqn += rho * canopy_->Fu();
    // Info << "canopyModelSource addRhoSup: " << fieldi << endl;
    // vectorField& Usource =  eqn.source();
    // auto Fu = canopy_->Fu();
    // Usource -= rho * mesh_.V() * Fu;
}

// Correcting all scalar equations
void Foam::fv::canopyModelSource::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{

    eqn += getScalarSource(eqn.psi().name());

    // Info << "canopyModelSource addSupS: "<< fieldi  << endl;
    // scalarField& source =  eqn.source();
    // auto Fi = getScalarSource(eqn.psi().name());
    // source -= Fi * mesh_.V();
        
}


void Foam::fv::canopyModelSource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
) 
{
    eqn += rho * getScalarSource(eqn.psi().name());

    // Info << "canopyModelSource addRhoSupC: " << eqn.psi().name() << endl;
    // scalarField& source =  eqn.source();
    // auto Fi = getScalarSource(eqn.psi().name());
    // source -= rho  * mesh_.V() * Fi;
}



void Foam::fv::canopyModelSource::correct(volVectorField& T)
{
    // Only update fvOption when the variable is T
    //if(T.name() == "T" && curTimeIndex_ != this->mesh().time().timeIndex())
    //{
      curTimeIndex_ = this->mesh().time().timeIndex();
      Info << "TreeModelSource Correcting.." << endl;
      canopy_->correctMomentumTransfer();
      canopy_->correctEnergyTransfer();
    //}
}

bool Foam::fv::canopyModelSource::read(const dictionary& dict)
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
