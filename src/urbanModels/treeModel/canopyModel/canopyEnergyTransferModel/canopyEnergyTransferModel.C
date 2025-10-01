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
    Base class to provide momentumTransport (include velocity, and associated turbulence
    parameters. 
    
SourceFiles
    canopyMomentumTransferModel.C
\*---------------------------------------------------------------------------*/

#include "canopyEnergyTransferModel.H"
#include "dimensionSets.H"
#include "dimensionedScalarFwd.H"
#include "scalar.H"


namespace Foam
{
template<class BaseCanopyModel>
void canopyEnergyTransferModel<BaseCanopyModel>::init()
{
    Info << "Empty init function" << endl;
}

template<class BaseCanopyModel>
canopyEnergyTransferModel<BaseCanopyModel>::canopyEnergyTransferModel
(
      const treeModel& tree
)
:
    BaseCanopyModel(tree),
    energyTransferDict_(BaseCanopyModel::dict().subDict("energyTransfer")),
    // radiation_(tree.mesh().lookupObjectRef<radiation::radiationModel>("radiationProperties")),
    // dom_(tree.mesh().lookupObjectRef<radiation::fvDOM>("radiationProperties")),
    FT_
    (
      IOobject
      (
        "canopy_FT",
        tree.mesh().time().constant(),
        tree.mesh(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
      ),
      tree.mesh(),
      dimensionedScalar(dimTemperature/dimTime, 0)
    ),
    Fq_
    (
      IOobject
      (
        "canopy_Fq",
        tree.mesh().time().constant(),
        tree.mesh(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
      ),
      tree.mesh(),
      dimensionedScalar(dimless/dimTime, 0)
    ),
    a_
    (
      IOobject
      (
        "canopy_a",
        tree.mesh().time().constant(),
        tree.mesh(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
      ),
      tree.mesh(),
      dimensionedScalar(dimless/dimLength, 0)
    ),
    e_
    (
      IOobject
      (
        "canopy_e",
        tree.mesh().time().constant(),
        tree.mesh(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
      ),
      tree.mesh(),
      dimensionedScalar(dimless/dimLength, 0)
    ),
    E_    
    (
      IOobject
      (
        "canopy_E",
        tree.mesh().time().constant(),
        tree.mesh(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
      ),
      tree.mesh(),
      dimensionedScalar(dimMass/dimLength/pow3(dimTime), 0)
    ),
    Tleaf_
    (
      IOobject
      (
        "canopy_Tleaf",
        tree.mesh().time().constant(),
        tree.mesh(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
      ),
      tree.mesh().lookupObjectRef<volScalarField>("T")
    ),
    Cleaf_
    (
        dimEnergy/dimMass/dimTemperature, 
        energyTransferDict_.get<double>("Cleaf")
    ),
    hleaf_
    (
        dimLength, 
        energyTransferDict_.get<double>("hleaf")
    ),
    rholeaf_
    (
        dimDensity, 
        energyTransferDict_.get<double>("rholeaf")
    ),
    rb_
    (
      IOobject
      (
        "canopy_rb",
        tree.mesh().time().constant(),
        tree.mesh(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
      ),
      tree.mesh(),
      dimensionedScalar(dimDensity, energyTransferDict_.get<double>("rb"))// copy the background temperature if not provided.
    ),
    rs_
    (
      IOobject
      (
        "canopy_rb",
        tree.mesh().time().constant(),
        tree.mesh(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
      ),
      tree.mesh(),
      dimensionedScalar(dimDensity, energyTransferDict_.get<double>("rs"))// copy the background temperature if not provided.
    )
{
    Info << "Testing instantiation of canopyEnergyTransferModel." << endl;
    // Info << "radiationType: " << radiation_.typeName_() << endl;
    // Info << "radiationType: " << radiation_.name() << endl;
    // Info << "radiationType: " << radiation_.type() << endl;

    init();
    correctEnergyTransfer();
};

template<class BaseCanopyModel>
void canopyEnergyTransferModel<BaseCanopyModel>::correctEnergyTransfer()
{
    Info << "Correcting energy transfer... " << endl;

    // Loading references to variables before correction iterations
    const fvMesh& mesh = BaseCanopyModel::tree().mesh();
    scalar deltaT = mesh.time().deltaT().value();
    // const volScalarField& rho = thermo_.rho()();
    // const volScalarField& Cp = thermo_.Cp()();
    // const volScalarField& T = thermo_.T();
    // const volScalarField& p = thermo_.p();
    scalar rho = 1.0;
    scalar Cp = 1005.0;
    scalar Lv = 2458.3;
    const volScalarField& T = mesh.lookupObjectRef<volScalarField>("T");
    const volScalarField& p = T.mesh().lookupObjectRef<volScalarField>("p");
    const volVectorField& U = BaseCanopyModel::U();
    const volScalarField& G = mesh.lookupObjectRef<volScalarField>("G");

    // Moisture related variables
    const volScalarField& q = T.mesh().lookupObjectRef<volScalarField>("H2O");

    // Setting coefficients relevant to the correction models
    // Stomatal resistance model
    scalar C_rb = 200;
    scalar rsmean = 242; // Maize from Turner & Begg, 1973
    scalar phi_rs = 40; // Meyers & Paw, 1987
    scalar rfcLeaf = 0.6; // leaf reflectance (arbitrary value)
    scalar C_ru = 0.5; // correct for Ru calculation (different direction)
                       //

    // Attempt for non-looping procedures
    // auto ldia = BaseCanopyModel::ldia();
    // auto lad = BaseCanopyModel::lad();
    // auto laCov = BaseCanopyModel::laCov();
    // dimensionedVector Usmall(dimVelocity, SMALL*vector(1,1,1));
    // dimensionedScalar Cpd(dimEnergy/dimTemperature, Cp);
    // dimensionedScalar rhod(dimDensity, rho);
    // dimensionedScalar Lvd(dimEnergy, Lv);
    // dimensionedScalar rsmeand(dimTime/dimLength, 242);
    // dimensionedScalar C_rbd(sqrt(dimTime), 200);
    // scalar C_ru = 0.5; // correct for Ru calculation (different direction)
    // dimensionedScalar rfcLeaf(dimVolume, 0.6); // leaf reflectance (arbitrary value)

    // // dimensionedScalar G(dimEnergy/dimArea, 200);
    // dimensionedScalar phi_rsd({1,0,-3,0,0}, 40);
    // 
    //   //- Aerodynamic resistance (Meyers & Paw, 1987)
    //   auto rbleaf = (C_rbd * sqrt(ldia / mag(U+Usmall)))();
    //   
    //   //- Convective heat transport from leaves to air
    //   auto FheConvd = (2 * lad * Cpd * (T - Tleaf_) / rbleaf)();
    //   //- Stomatal resistance (Meyers & Paw, 1987)
    //   auto rsleaf = rsmeand * (1 + phi_rsd / (G + dimensionedScalar({1,0,-3,0,0}, SMALL)));
    //   // Info << "test 1" << endl;
    //   
    //   //- Bolton 1980 for saturation vapor pressure
    //   //- Convert temperature to degree celcius
    //   auto pTcelleaf = Tleaf_ - dimensionedScalar(dimTemperature, 273.15);
    //   auto Tcelleaf = pTcelleaf.ref();
    //   
    //   auto esatleaf
    //   (
    //     (
    //       dimensionedScalar(dimPressure/dimDensity, 100) * 6.112 * exp
    //       (
    //         Tcelleaf * 17.67 / (Tcelleaf + dimensionedScalar(dimTemperature , 243.5))
    //       )
    //     )()
    //   );
    //   auto qsatleaf
    //   (
    //     (
    //       (0.622 * esatleaf) / (p - (0.378 * esatleaf))
    //     )()
    //   );
    //   // Info << "test 4" << endl;
    //   // Info << qsatleaf.dimensions() << q.dimensions() << endl;
    //   // Info << rbleaf.dimensions() << rsleaf.dimensions() << endl;
    //   auto Fqleaf 
    //   (
    //       (2 * lad * (qsatleaf - q)  / (rbleaf + rsleaf))()
    //   ); // Placeholder 

    //   // Info << "test 5" << endl;
    //   auto FheLeConvd = (Lvd * Fqleaf)();
    //   // Info << average(FheLeConvd) << endl;
    //   // Info << "test 6" << endl;

    //  //- Radiative Transfer
    //    // Need a conversion function for total leaf area and absorption
    //    // coefficient
    //    // This version is based on a mapping of leaf coverage ratio to 
    //    // absorption coefficient

    //    auto aleaf = ((pos(lad) * 0.5) / dimensionedScalar(dimLength, 1))();
    //    auto eleaf = ((pos(lad) * 0.5) / dimensionedScalar(dimLength, 1))();
    //   Info << "aleaf: " << average(aleaf) << endl;
    //   
    //    auto Eleaf = (eleaf * constant::physicoChemical::sigma * pow4(Tleaf_) * pos(lad) )();
    //    auto Ruleaf = (C_ru * rfcLeaf * laCov * G)() ;

    //    a_ = aleaf;
    //    e_ = eleaf;
    //    E_ = Eleaf;

    //    Info << "Eleaf: " << Eleaf.dimensions() << endl;
    //    Info << "Ruleaf: " << Ruleaf.dimensions() << endl;
    //    Info << "FheConvd: " << FheConvd.dimensions() << endl;
    //    // auto dleaf = (Ruleaf - (FheConvd + FheLeConvd + Eleaf * dimensionedScalar(dimVolume, 1.0))) / (Cleaf_ * rholeaf_ * lad * hleaf_);
    //    auto m_leaf = (Cleaf_ * rholeaf_ * lad * hleaf_)();
    //    Info << "m_leaf: " << m_leaf.dimensions() << endl;
    //    forAll(m_leaf, i)
    //    {
    //       auto mi = m_leaf[i];
    //       if(mi <= 0) 
    //       {
    //           m_leaf.primitiveFieldRef()[i] = GREAT;
    //       }
    //    }
    //    Info << "m_leaf: " << gMin(m_leaf) << endl;
    //    // Units for dleaf is way off, need to harmonize this 
    //    auto dleaf = (Ruleaf - (FheConvd + FheLeConvd)) / m_leaf;

    //   // Correcting Energy Related Terms
    //   // FT_[i] = (-1) * (FheConv + FheLeConv) / (CpCell * rhoCell);
    //   FT_ = (-1) * (FheConvd + FheLeConvd ) / Cpd;
    //   Fq_ = Fqleaf;
    //   Tleaf_.primitiveFieldRef() = (Tleaf_ + dleaf * deltaT)->primitiveField();

    // Pout << "Starting iter on " << mesh.V().size() << " cells"  << endl;
    forAll(mesh.V(), i)
    {
        scalar ladCell = BaseCanopyModel::lad()[i];
        scalar ldiaCell = BaseCanopyModel::ldia()[i];
        scalar laleafCell = BaseCanopyModel::la()[i];
        scalar laCovCell = BaseCanopyModel::laCov()[i];
        // Info << i << ", " << laleafCell << endl;

        if (laleafCell <= 0){continue; }
        scalar TleafCell = Tleaf_[i];
        scalar GCell = G[i];
        // scalar GCell = 200.0;
        // Loading cell values for previous loaded variables
        scalar rhoCell = rho;
        scalar CpCell = Cp;
        vector UCell = U[i];
        scalar TCell = T[i];
        scalar pCell = p[i];
        scalar qCell = q[i];
        scalar rholeafCell = rholeaf_.value();
        scalar CleafCell = Cleaf_.value();
        scalar hleafCell = hleaf_.value();


        // Correcting dependent variables

          //- Aerodynamic resistance (Meyers & Paw, 1987)
          scalar rbleafCell = C_rb * 
              sqrt
              (
                ldiaCell / mag(UCell+SMALL*vector(1,1,1))
              );
          
          //- Convective heat transport from leaves to air
          scalar FheConv = 2 * ladCell * CpCell * rhoCell * (TleafCell - TCell) / rbleafCell;
          //- Stomatal resistance (Meyers & Paw, 1987)
          scalar rsleafCell = rsmean * (1 + phi_rs / (GCell+1));
          
          //- Bolton 1980 for saturation vapor pressure
          //- Convert temperature to degree celcius
          scalar TcelleafCell = TleafCell - 273.15;
          
          scalar esatleafCell
          (
            100 * 6.112 * exp
            (
              TcelleafCell * 17.67 / (TcelleafCell + 243.5)
            )
          );
          scalar qsatleafCell
          (
            (0.622 * esatleafCell) / (pCell - (0.378 * esatleafCell))
          );
          scalar Fqleaf 
          (
              2 * ladCell * (qsatleafCell - qCell)  / (rbleafCell + rsleafCell)

          ); // Placeholder 
          Fqleaf *= pos(Fqleaf);
          // Info << ladCell << ", " 
          //      << qsatleafCell << ", "
          //      << qCell << ", "
          //      << rbleafCell << ", "
          //      << rsleafCell << ", "
          //      << endl;
          scalar FheLeConv = rhoCell * Lv * Fqleaf;

        //- Radiative Transfer
          // Need a conversion function for total leaf area and absorption
          // coefficient
          // This version is based on a mapping of leaf coverage ratio to 
          // absorption coefficient
          // scalar aleafCell = (-1) * log(1 - min(laCovCell, 0.9999));
          // scalar eleafCell = aleafCell * 0.8;
          // Not using laLit as a is irrespective of leaves being lit
          scalar aleafCell = ladCell * 0.1;
          scalar eleafCell = ladCell * 0.1;

          scalar EleafCell = eleafCell * constant::physicoChemical::sigma.value() * pow4(TleafCell);

          //- Radiative heat input from incidental radiative fluxes
          //- Note this assumes isotropic leaf absorption
          //- Since the cell leaves are not spherical, C_ru is applied to
          // reduce the total amount of incident radiation. 
          //- rfcLeaf is applied to account for the reflectance of leaves, 
          //  however, it should be noted that reflectance should be correctly
          //  accounted for in the scattering model.
          scalar Ruleaf = C_ru * rfcLeaf * laCovCell * GCell ;

        scalar dleaf =  (Ruleaf - (FheConv + FheLeConv + EleafCell*0)) /  (CleafCell * rholeafCell * ladCell * hleafCell);
          // Info << Ruleaf << ", "
          //      << FheConv << ", " << TleafCell << ", " << TCell << ", " << rbleafCell << ", " << ladCell << ", "
          //      << FheLeConv << ", "
          //      << EleafCell << ", "
          //      << endl;

        // Update field values
        a_[i] = aleafCell;
        e_[i] = eleafCell;
        E_[i] = EleafCell;

        rb_[i] = rbleafCell;
        rs_[i] = rsleafCell;
        
        // Correcting Energy Related Terms
        scalar m = 1.0;
        FT_[i] = (-1 * m ) * (FheConv + FheLeConv) / (CpCell * rhoCell);
        // FT_[i] = (-1) * (FheConv ) / (CpCell * rhoCell);
        Fq_[i] = m * Fqleaf;
        Tleaf_[i] = TleafCell + m * dleaf * deltaT;
    }
}

}
