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
    Foam::canopyMomentumTransferModel

Description
    Base class to provide thermo-related transfer model (include temperature, and)

    Leaf properties names are referenced from here (partially): 

    https://www.researchgate.net/figure/Thermal-properties-of-seven-species-of-plant-leaves_tbl1_215452202#:~:text=The%20specific%20heat%20capacity%20of%20fresh%20leaves%20varies%20from%201255,leaves%20is%20necessary.%20...

SourceFiles
    canopyMomentumTransferModel.C
\*---------------------------------------------------------------------------*/

#ifndef canopyEnergyTransferModel_H
#define canopyEnergyTransferModel_H

#include "treeModel.H"
// #include "fluidAtmThermo.H"
#include "radiationModel.H"
#include "fvDOM.H"
#include "cellSetHashTables.H"

namespace Foam
{

template<class BaseCanopyModel>
class canopyEnergyTransferModel
  :
      public BaseCanopyModel
{

protected:
    
    
        //- Dictionary for momentumTransfer
        dictionary energyTransferDict_;

        //- Reference to Radiation model
        // radiation::radiationModel& radiation_;
        // radiation::fvDOM& dom_;


        volScalarField FT_;
        volScalarField Fq_;

        // Radiation related
        volScalarField a_;
        volScalarField e_;
        volScalarField E_;

        //- Properties of leaves of cell
        
          // Temperature of leaves
          volScalarField  Tleaf_;
          
          // Heat capacity of leaves
          dimensionedScalar  Cleaf_;

          // Leaf thickness
          dimensionedScalar  hleaf_;

          // Leaf density
          dimensionedScalar  rholeaf_;

          // Leaf Aerodynamic Resistance
          volScalarField  rb_;

          // Leaf Stomatal Resistance
          volScalarField  rs_;

          

private: 
          void init();

public:
      
    canopyEnergyTransferModel
    (
      const treeModel& tree
    );

    virtual ~canopyEnergyTransferModel(){};

    //- Member Functions
    virtual volScalarField& FT(){return FT_;};
    virtual const volScalarField& FT() const{return FT_;};
    virtual volScalarField& Fq(){return Fq_;};
    virtual const volScalarField& Fq() const{return Fq_;};
    virtual volScalarField& Tleaf() {return Tleaf_;};
    virtual const volScalarField& Tleaf() const {return Tleaf_;};

    //- Radiation related 
    virtual volScalarField& a(){return a_;};
    virtual const volScalarField& a() const{return a_;};
    virtual volScalarField& e(){return e_;};
    virtual const volScalarField& e() const{return e_;};
    virtual volScalarField& E(){return E_;};
    virtual const volScalarField& E() const{return E_;};

    virtual void correctEnergyTransfer();
};

}

#ifdef NoRepository
  #include "canopyEnergyTransferModel.C"
#endif

#endif

