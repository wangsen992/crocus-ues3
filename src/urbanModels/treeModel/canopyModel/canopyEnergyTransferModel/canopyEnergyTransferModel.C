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


namespace Foam
{
template<class BaseCanopyModel>
void canopyEnergyTransferModel<BaseCanopyModel>::init()
{
    Info << "Initiating Tleaf and Cleaf" << endl;

    dimensionedScalar CleafDefault
    (
        dimEnergy/dimMass/dimTemperature, 
        energyTransferDict_.lookup<double>("Cleaf")
    );
    dimensionedScalar hleafDefault
    (
        dimLength, 
        energyTransferDict_.lookup<double>("hleaf")
    );
    dimensionedScalar rholeafDefault
    (
        dimDensity, 
        energyTransferDict_.lookup<double>("rholeaf")
    );
    forAllConstIter(labelHashSet, BaseCanopyModel::canopyCells(), iter)
    {
        Fhe_.set
        (
          iter.key(),
          dimensionedScalar(dimEnergy/dimTime, 0)
        );
        Fq_.set
        (
          iter.key(),
          dimensionedScalar(dimless, 0)
        );
        a_.set
        (
          iter.key(),
          dimensionedScalar(dimless/dimLength, 0)
        );
        e_.set
        (
          iter.key(),
          dimensionedScalar(dimless/dimLength, 0)
        );
        E_.set
        (
          iter.key(),
          dimensionedScalar(dimMass/dimLength/pow3(dimTime), 0)
        );
        Tleaf_.set
        (
          iter.key(),
          dimensionedScalar(dimTemperature, thermo_.T()[iter.key()])
        );
        Cleaf_.set
        (
          iter.key(),
          CleafDefault
        );
        hleaf_.set
        (
          iter.key(),
          hleafDefault
        );
        rholeaf_.set
        (
          iter.key(),
          rholeafDefault
        );
        rb_.set
        (
          iter.key(),
          dimensionedScalar(dimless/dimVelocity, 2)
        );
        rs_.set
        (
          iter.key(),
          dimensionedScalar(dimless/dimVelocity, 2)
        );
    }
}

template<class BaseCanopyModel>
canopyEnergyTransferModel<BaseCanopyModel>::canopyEnergyTransferModel
(
      const treeModel& tree
)
:
    BaseCanopyModel(tree),
    energyTransferDict_(BaseCanopyModel::dict().subDict("energyTransfer")),
    thermo_(tree.mesh().lookupObjectRef<fluidAtmThermo>(fluidAtmThermo::dictName)),
    radiation_(tree.mesh().lookupObjectRef<radiationModel>("radiationProperties")),
    dom_(tree.mesh().lookupObjectRef<radiationModels::fvDOM>("radiationProperties")),
    Fhe_(BaseCanopyModel::canopyCells().size()),
    Fq_(BaseCanopyModel::canopyCells().size()),
    a_(BaseCanopyModel::canopyCells().size()),
    e_(BaseCanopyModel::canopyCells().size()),
    E_(BaseCanopyModel::canopyCells().size()),
    Tleaf_(BaseCanopyModel::canopyCells().size()),
    Cleaf_(BaseCanopyModel::canopyCells().size()),
    hleaf_(BaseCanopyModel::canopyCells().size()),
    rholeaf_(BaseCanopyModel::canopyCells().size()),
    rb_(BaseCanopyModel::canopyCells().size()),
    rs_(BaseCanopyModel::canopyCells().size())
{
    Info << "Testing instantiation of canopyEnergyTransferModel." << endl;
    Info << "thermoType: " << thermo_.thermoName() << endl;
    Info << "radiationType: " << radiation_.typeName_() << endl;
    Info << "radiationType: " << radiation_.name() << endl;
    Info << "radiationType: " << radiation_.type() << endl;

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
    const volScalarField& rho = thermo_.rho()();
    const volScalarField& Cp = thermo_.Cp()();
    const volScalarField& T = thermo_.T();
    const volScalarField& p = thermo_.p();
    const volVectorField& U = BaseCanopyModel::U();
    const volScalarField& G = dom_.G();

    // Moisture related variables
    const volScalarField& q = thermo_.composition().Y("H2O");
    label qIndex = thermo_.composition().index(q);
    scalar Lv = thermo_.composition().Hf(qIndex);

    // Setting coefficients relevant to the correction models
    // Stomatal resistance model
    scalar C_rb = 200;
    scalar rsmean = 242; // Maize from Turner & Begg, 1973
    scalar phi_rs = 40; // Meyers & Paw, 1987
    scalar rfcLeaf = 0.6; // leaf reflectance (arbitrary value)
    scalar C_ru = 0.5; // correct for Ru calculation (different direction)

    Pout << "Starting iter on " << BaseCanopyModel::canopyCells().size() << " cells"  << endl;
    forAllConstIter(labelHashSet, BaseCanopyModel::canopyCells(), iter)
    {
        // Loading cell values for previous loaded variables
        vector UCell = U[iter.key()];
        scalar rhoCell = rho[iter.key()];
        scalar CpCell = Cp[iter.key()];
        scalar TCell = T[iter.key()];
        scalar pCell = p[iter.key()];
        scalar qCell = q[iter.key()];
        scalar TleafCell = Tleaf_[iter.key()].value();
        scalar CleafCell = Cleaf_[iter.key()].value();
        scalar hleafCell = hleaf_[iter.key()].value();
        scalar rholeafCell = rholeaf_[iter.key()].value();
        scalar GCell = G[iter.key()];
        scalar laleafCell = BaseCanopyModel::la()[iter.key()].value();
        scalar laCovCell = BaseCanopyModel::laCov()[iter.key()].value();
        scalar ldiaCell = BaseCanopyModel::ldia()[iter.key()].value();

        // Correcting dependent variables

          //- Aerodynamic resistance (Meyers & Paw, 1987)
          scalar rbleafCell = C_rb * 
              sqrt
              (
                ldiaCell / mag(UCell)
              );
          
          //- Convective heat transport from leaves to air
          scalar FheConv = 2 * laleafCell * CpCell * rhoCell * (TleafCell - TCell) / rbleafCell;
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
            2 * laleafCell * (qsatleafCell - qCell)  / (rbleafCell + rsleafCell)
          ); // Placeholder 
          scalar FheLeConv = rhoCell * Lv *  Fqleaf;

        //- Radiative Transfer
          // Need a conversion function for total leaf area and absorption
          // coefficient
          // This version is based on a mapping of leaf coverage ratio to 
          // absorption coefficient
          scalar aleafCell = (-1) * log(1 - min(laCovCell, 0.9999));
          scalar eleafCell = aleafCell * 0.8;
          // Not using laLit as a is irrespective of leaves being lit
          // scalar aleafCell = BaseCanopyModel::laLit()[iter.key()].value() / 0.15;
          // scalar eleafCell = BaseCanopyModel::laLit()[iter.key()].value() / 0.15;

          scalar EleafCell = eleafCell * constant::physicoChemical::sigma.value() * pow4(TleafCell);

          //- Radiative heat input from incidental radiative fluxes
          //- Note this assumes isotropic leaf absorption
          //- Since the cell leaves are not spherical, C_ru is applied to
          // reduce the total amount of incident radiation. 
          //- rfcLeaf is applied to account for the reflectance of leaves, 
          //  however, it should be noted that reflectance should be correctly
          //  accounted for in the scattering model.
          scalar Ruleaf = C_ru * rfcLeaf * laCovCell * GCell ;

        scalar dleaf =  (Ruleaf - (FheConv + FheLeConv + EleafCell)) /  (CleafCell * rholeafCell * laleafCell * hleafCell);

        // Pout << Ruleaf << " - " << "( " << FheConv << " + "
        //      << FheLeConv << " + " << EleafCell << " ) " 
        //      << " = " << dleaf << endl;

        // Update field values
        a_[iter.key()].value() = aleafCell;
        e_[iter.key()].value() = eleafCell;
        E_[iter.key()].value() = EleafCell;

        rb_[iter.key()].value() = rbleafCell;
        rs_[iter.key()].value() = rsleafCell;
        
        // Correcting Energy Related Terms
        Fhe_[iter.key()].value() = FheConv;
        Fq_[iter.key()].value() = Fqleaf;
        Tleaf_[iter.key()].value() = TleafCell + dleaf * deltaT;
    }
}

}
