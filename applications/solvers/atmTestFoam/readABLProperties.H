    IOdictionary ABLProperties
    (
        IOobject
        (
            "ABLProperties",
            runTime.time().constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // PROPERTIES CONCERNING THE SOURCE TERMS.

       // Specify the type of source to use for momentum and temperature.  The
       // possible types are "given" and "computed".  
       // - The "given" type means that the source values are directly given
       //   and the momentum and temperature fields will react accordingly.  
       // - The "computed" type means that the mean velocity and temperature
       //   are given and the source terms that maintain them are computed. 
       word momentumSourceType(ABLProperties.lookup("momentumSourceType"));
       word temperatureSourceType(ABLProperties.lookup("temperatureSourceType"));

       // If giving the velocity and computing the sources, specify how the velocity
       // is given.  "component" means you enter the x, y, anc z components.
       // "speedAndDirection" means that you enter the horizontal wind speed, horizontal
       // direction, and vertical component.
       word velocityInputType(ABLProperties.lookup("velocityInputType"));

       // Read in the heights at which the sources are given.
       List<scalar> sourceHeightsMomentumSpecified(ABLProperties.lookup("Z"));
       List<scalar> sourceHeightsTemperatureSpecified(ABLProperties.lookup("Z"));

       // Read in LS2D Terms
       List<List<scalar> > source_dudt_advec(ABLProperties.lookup("dudt_advec"));
       List<List<scalar> > source_dvdt_advec(ABLProperties.lookup("dvdt_advec"));
       List<List<scalar> > source_dthldt_advec(ABLProperties.lookup("dthldt_advec"));
       List<List<scalar> > source_dqtdt_advec(ABLProperties.lookup("dqtdt_advec"));
       List<List<scalar> > source_u(ABLProperties.lookup("u"));
       List<List<scalar> > source_v(ABLProperties.lookup("v"));
       List<List<scalar> > source_w(ABLProperties.lookup("w"));
       List<List<scalar> > source_ug(ABLProperties.lookup("ug"));
       List<List<scalar> > source_vg(ABLProperties.lookup("vg"));
       List<List<scalar> > source_thl(ABLProperties.lookup("thl"));
       List<List<scalar> > source_qt(ABLProperties.lookup("qt"));

       // Read in the table of momentum source vs. time and height.
       List<List<scalar> > sourceTableMomentumX(ABLProperties.lookup("dudt_advec"));
       List<List<scalar> > sourceTableMomentumY(ABLProperties.lookup("dvdt_advec"));
       List<List<scalar> > sourceTableMomentumZ(ABLProperties.lookup("w"));

       // Read in the table of temperature source vs. time and height.
       List<List<scalar> > sourceTableTemperature(ABLProperties.lookup("dthldt_advec"));


       // Check sizes of source tables.
       label nSourceMomentumHeights = sourceHeightsMomentumSpecified.size();
       label nSourceTemperatureHeights = sourceHeightsTemperatureSpecified.size();

       label nSource_dudt_advec_Times = source_dudt_advec.size();
       label nSource_dvdt_advec_Times = source_dvdt_advec.size();
       label nSource_dthldt_advec_Times = source_dthldt_advec.size();
       label nSource_dqtdt_advec_Times = source_dqtdt_advec.size();
       label nSource_u_Times = source_u.size();
       label nSource_v_Times = source_v.size();
       label nSource_w_Times = source_w.size();
       label nSource_ug_Times = source_ug.size();
       label nSource_vg_Times = source_vg.size();
       label nSource_thl_Times = source_thl.size();
       label nSource_qt_Times = source_qt.size();

       label nSourceMomentumXTimes = sourceTableMomentumX.size();
       label nSourceMomentumYTimes = sourceTableMomentumY.size();
       label nSourceMomentumZTimes = sourceTableMomentumZ.size();
       label nSourceTemperatureTimes = sourceTableTemperature.size();

       forAll(source_dudt_advec,i)
       {
           if (source_dudt_advec[i].size()-1 != nSourceMomentumHeights)
           {
               FatalErrorIn
               (
                   "Number of source_dudt_advec heights does not match sourceHeightsMomentum"
               )   << abort(FatalError);
           }
       }
       forAll(source_dvdt_advec,i)
       {
           if (source_dvdt_advec[i].size()-1 != nSourceMomentumHeights)
           {
               FatalErrorIn
               (
                   "Number of source_dvdt_advec heights does not match sourceHeightsMomentum"
               )   << abort(FatalError);
           }
       }
       forAll(source_dthldt_advec,i)
       {
           if (source_dthldt_advec[i].size()-1 != nSourceMomentumHeights)
           {
               FatalErrorIn
               (
                   "Number of source_dthldt_advec heights does not match sourceHeightsMomentum"
               )   << abort(FatalError);
           }
       }
       forAll(source_dqtdt_advec,i)
       {
           if (source_dqtdt_advec[i].size()-1 != nSourceMomentumHeights)
           {
               FatalErrorIn
               (
                   "Number of source_dqtdt_advec heights does not match sourceHeightsMomentum"
               )   << abort(FatalError);
           }
       }
       forAll(source_u,i)
       {
           if (source_u[i].size()-1 != nSourceMomentumHeights)
           {
               FatalErrorIn
               (
                   "Number of source_u heights does not match sourceHeightsMomentum"
               )   << abort(FatalError);
           }
       }
       forAll(source_v,i)
       {
           if (source_v[i].size()-1 != nSourceMomentumHeights)
           {
               FatalErrorIn
               (
                   "Number of source_v heights does not match sourceHeightsMomentum"
               )   << abort(FatalError);
           }
       }
       forAll(source_w,i)
       {
           if (source_w[i].size()-1 != nSourceMomentumHeights)
           {
               FatalErrorIn
               (
                   "Number of source_w heights does not match sourceHeightsMomentum"
               )   << abort(FatalError);
           }
       }
       forAll(source_ug,i)
       {
           if (source_ug[i].size()-1 != nSourceMomentumHeights)
           {
               FatalErrorIn
               (
                   "Number of source_ug heights does not match sourceHeightsMomentum"
               )   << abort(FatalError);
           }
       }
       forAll(source_vg,i)
       {
           if (source_vg[i].size()-1 != nSourceMomentumHeights)
           {
               FatalErrorIn
               (
                   "Number of source_vg heights does not match sourceHeightsMomentum"
               )   << abort(FatalError);
           }
       }
       forAll(source_thl,i)
       {
           if (source_thl[i].size()-1 != nSourceMomentumHeights)
           {
               FatalErrorIn
               (
                   "Number of source_vg heights does not match sourceHeightsMomentum"
               )   << abort(FatalError);
           }
       }
       forAll(source_qt,i)
       {
           if (source_qt[i].size()-1 != nSourceMomentumHeights)
           {
               FatalErrorIn
               (
                   "Number of source_vg heights does not match sourceHeightsMomentum"
               )   << abort(FatalError);
           }
       }
       forAll(sourceTableMomentumX,i)
       {
           if (sourceTableMomentumX[i].size()-1 != nSourceMomentumHeights)
           {
               FatalErrorIn
               (
                   "Number of sourceTableMomentumX heights does not match sourceHeightsMomentum"
               )   << abort(FatalError);
           }
       }
       forAll(sourceTableMomentumY,i)
       {
           if (sourceTableMomentumY[i].size()-1 != nSourceMomentumHeights)
           {
               FatalErrorIn
               (
                   "Number of sourceTableMomentumY heights does not match sourceHeightsMomentum"
               )   << abort(FatalError);
           }
       }
       forAll(sourceTableMomentumZ,i)
       {
           if (sourceTableMomentumZ[i].size()-1 != nSourceMomentumHeights)
           {
               FatalErrorIn
               (
                   "Number of sourceTableMomentumZ heights does not match sourceHeightsMomentum"
               )   << abort(FatalError);
           }
       }
       forAll(sourceTableTemperature,i)
       {
           if (sourceTableTemperature[i].size()-1 != nSourceTemperatureHeights)
           {
               FatalErrorIn
               (
                   "Number of sourceTableTemperature heights does not match sourceHeightsTemperature"
               )   << abort(FatalError);
           }
       }


       // Break the source tables into interpolation tables.
       List<scalar> source_dudt_advecTimesSpecified(nSource_dudt_advec_Times,0.0);
       List<scalar> source_dvdt_advecTimesSpecified(nSource_dvdt_advec_Times,0.0);
       List<scalar> source_dthldt_advecTimesSpecified(nSource_dthldt_advec_Times,0.0);
       List<scalar> source_dqtdt_advecTimesSpecified(nSource_dqtdt_advec_Times,0.0);
       List<scalar> source_uTimesSpecified(nSource_u_Times,0.0);
       List<scalar> source_vTimesSpecified(nSource_v_Times,0.0);
       List<scalar> source_wTimesSpecified(nSource_w_Times,0.0);
       List<scalar> source_ugTimesSpecified(nSource_ug_Times,0.0);
       List<scalar> source_vgTimesSpecified(nSource_vg_Times,0.0);
       List<scalar> source_thlTimesSpecified(nSource_thl_Times,0.0);
       List<scalar> source_qtTimesSpecified(nSource_qt_Times,0.0);

       List<scalar> sourceMomentumXTimesSpecified(nSourceMomentumXTimes,0.0);
       List<scalar> sourceMomentumYTimesSpecified(nSourceMomentumYTimes,0.0);
       List<scalar> sourceMomentumZTimesSpecified(nSourceMomentumZTimes,0.0);
       List<scalar> sourceTemperatureTimesSpecified(nSourceTemperatureTimes,0.0);


       List<List<scalar> > source_dudt_advecSpecified(nSource_dudt_advec_Times,List<scalar>(nSourceMomentumHeights,0.0));
       List<List<scalar> > source_dvdt_advecSpecified(nSource_dvdt_advec_Times,List<scalar>(nSourceMomentumHeights,0.0));
       List<List<scalar> > source_dthldt_advecSpecified(nSource_dthldt_advec_Times,List<scalar>(nSourceMomentumHeights,0.0));
       List<List<scalar> > source_dqtdt_advecSpecified(nSource_dqtdt_advec_Times,List<scalar>(nSourceMomentumHeights,0.0));
       List<List<scalar> > source_uSpecified(nSource_u_Times,List<scalar>(nSourceMomentumHeights,0.0));
       List<List<scalar> > source_vSpecified(nSource_v_Times,List<scalar>(nSourceMomentumHeights,0.0));
       List<List<scalar> > source_wSpecified(nSource_w_Times,List<scalar>(nSourceMomentumHeights,0.0));
       List<List<scalar> > source_ugSpecified(nSource_ug_Times,List<scalar>(nSourceMomentumHeights,0.0));
       List<List<scalar> > source_vgSpecified(nSource_vg_Times,List<scalar>(nSourceMomentumHeights,0.0));
       List<List<scalar> > source_thlSpecified(nSource_thl_Times,List<scalar>(nSourceMomentumHeights,0.0));
       List<List<scalar> > source_qtSpecified(nSource_qt_Times,List<scalar>(nSourceMomentumHeights,0.0));

       List<List<scalar> > sourceMomentumXSpecified(nSourceMomentumXTimes,List<scalar>(nSourceMomentumHeights,0.0));
       List<List<scalar> > sourceMomentumYSpecified(nSourceMomentumYTimes,List<scalar>(nSourceMomentumHeights,0.0));
       List<List<scalar> > sourceMomentumZSpecified(nSourceMomentumZTimes,List<scalar>(nSourceMomentumHeights,0.0));
       List<List<scalar> > sourceTemperatureSpecified(nSourceTemperatureTimes,List<scalar>(nSourceTemperatureHeights,0.0));

       for(int i = 0; i < nSource_dudt_advec_Times; i++)
       {
           source_dudt_advecTimesSpecified[i] = source_dudt_advec[i][0];
           for(int j = 0; j < nSourceMomentumHeights; j++)
           {
               source_dudt_advecSpecified[i][j] = source_dudt_advec[i][j+1];
           }
       }
       for(int i = 0; i < nSource_dvdt_advec_Times; i++)
       {
           source_dvdt_advecTimesSpecified[i] = source_dvdt_advec[i][0];
           for(int j = 0; j < nSourceMomentumHeights; j++)
           {
               source_dvdt_advecSpecified[i][j] = source_dvdt_advec[i][j+1];
           }
       }
       for(int i = 0; i < nSource_dthldt_advec_Times; i++)
       {
           source_dthldt_advecTimesSpecified[i] = source_dthldt_advec[i][0];
           for(int j = 0; j < nSourceMomentumHeights; j++)
           {
               source_dthldt_advecSpecified[i][j] = source_dthldt_advec[i][j+1];
           }
       }
       for(int i = 0; i < nSource_dqtdt_advec_Times; i++)
       {
           source_dqtdt_advecTimesSpecified[i] = source_dqtdt_advec[i][0];
           for(int j = 0; j < nSourceMomentumHeights; j++)
           {
               source_dqtdt_advecSpecified[i][j] = source_dqtdt_advec[i][j+1];
           }
       }
       for(int i = 0; i < nSource_u_Times; i++)
       {
           source_uTimesSpecified[i] = source_u[i][0];
           for(int j = 0; j < nSourceMomentumHeights; j++)
           {
               source_uSpecified[i][j] = source_u[i][j+1];
           }
       }
       for(int i = 0; i < nSource_v_Times; i++)
       {
           source_vTimesSpecified[i] = source_v[i][0];
           for(int j = 0; j < nSourceMomentumHeights; j++)
           {
               source_vSpecified[i][j] = source_v[i][j+1];
           }
       }
       for(int i = 0; i < nSource_w_Times; i++)
       {
           source_wTimesSpecified[i] = source_w[i][0];
           for(int j = 0; j < nSourceMomentumHeights; j++)
           {
               source_wSpecified[i][j] = source_w[i][j+1];
           }
       }
       for(int i = 0; i < nSource_ug_Times; i++)
       {
           source_ugTimesSpecified[i] = source_ug[i][0];
           for(int j = 0; j < nSourceMomentumHeights; j++)
           {
               source_ugSpecified[i][j] = source_ug[i][j+1];
           }
       }
       for(int i = 0; i < nSource_vg_Times; i++)
       {
           source_vgTimesSpecified[i] = source_vg[i][0];
           for(int j = 0; j < nSourceMomentumHeights; j++)
           {
               source_vgSpecified[i][j] = source_vg[i][j+1];
           }
       }
       for(int i = 0; i < nSource_thl_Times; i++)
       {
           source_thlTimesSpecified[i] = source_thl[i][0];
           for(int j = 0; j < nSourceMomentumHeights; j++)
           {
               source_thlSpecified[i][j] = source_thl[i][j+1];
           }
       }
       for(int i = 0; i < nSource_qt_Times; i++)
       {
           source_qtTimesSpecified[i] = source_qt[i][0];
           for(int j = 0; j < nSourceMomentumHeights; j++)
           {
               source_qtSpecified[i][j] = source_qt[i][j+1];
           }
       }

       for(int i = 0; i < nSourceMomentumXTimes; i++)
       {
           sourceMomentumXTimesSpecified[i] = sourceTableMomentumX[i][0];
           for(int j = 0; j < nSourceMomentumHeights; j++)
           {
               sourceMomentumXSpecified[i][j] = sourceTableMomentumX[i][j+1];
           }
       }

       for(int i = 0; i < nSourceMomentumYTimes; i++)
       {
           sourceMomentumYTimesSpecified[i] = sourceTableMomentumY[i][0];
           for(int j = 0; j < nSourceMomentumHeights; j++)
           {
               sourceMomentumYSpecified[i][j] = sourceTableMomentumY[i][j+1];
           }
       }

       for(int i = 0; i < nSourceMomentumZTimes; i++)
       {
           sourceMomentumZTimesSpecified[i] = sourceTableMomentumZ[i][0];
           for(int j = 0; j < nSourceMomentumHeights; j++)
           {
               sourceMomentumZSpecified[i][j] = sourceTableMomentumZ[i][j+1];
           }
       }

       for(int i = 0; i < nSourceTemperatureTimes; i++)
       {
           sourceTemperatureTimesSpecified[i] = sourceTableTemperature[i][0];
           for(int j = 0; j < nSourceTemperatureHeights; j++)
           {
               sourceTemperatureSpecified[i][j] = sourceTableTemperature[i][j+1];
           }
       }


       // If the desired mean wind or temperature is given at only one height, then revert to
       // the old way of specifying the source term.  Find the two grid levels that bracket
       // the given level and interpolate in between them to compute the source term.  So,
       // here, find those two levels.
       label hLevelsWind1I = 0;
       label hLevelsWind2I = 0;
       scalar hLevelsWind1 = 0;
       scalar hLevelsWind2 = 0;
       if ((momentumSourceType == "computed") && (nSourceMomentumHeights == 1))
       {
           #include "findWindHeight.H"
       }

       label hLevelsTemp1I = 0;
       label hLevelsTemp2I = 0;
       scalar hLevelsTemp1 = 0;
       scalar hLevelsTemp2 = 0;
       if ((temperatureSourceType == "computed") && (nSourceTemperatureHeights == 1))
       {
           #include "findTemperatureHeight.H"
       }


       // Relaxation factor on the source term application.
       scalar alphaMomentum(ABLProperties.lookupOrDefault<scalar>("alphaMomentum",1.0));
       scalar alphaTemperature(ABLProperties.lookupOrDefault<scalar>("alphaTemperature",1.0));


       // Large scale nudging relaxation value
       scalar tau(ABLProperties.lookupOrDefault<scalar>("tau", 1e-4));


    // PROPERTIES CONCERNING CORIOLIS FORCES

       // Planetary rotation period (hours)
       scalar planetaryRotationPeriod(readScalar(ABLProperties.lookup("planetaryRotationPeriod")));

       // Latitude on the planetary body (degrees)
       scalar latitude(readScalar(ABLProperties.lookup("latitude")));

       // Compute the planetar rotation vector
       vector Omega_;
       Omega_.x() = 0.0;
       Omega_.y() = ((2.0 * Foam::constant::mathematical::pi) / (max(1.0E-5,planetaryRotationPeriod)*3600.0)) * Foam::cos(latitude*Foam::constant::mathematical::pi/180.0);
       Omega_.z() = ((2.0 * Foam::constant::mathematical::pi) / (max(1.0E-5,planetaryRotationPeriod)*3600.0)) * Foam::sin(latitude*Foam::constant::mathematical::pi/180.0);
       uniformDimensionedVectorField Omega
       (
           IOobject
           (
               "Omega",
               runTime.constant(),
               mesh,
               IOobject::NO_READ,
               IOobject::NO_WRITE
           ),
           dimensionedVector("Omega",dimensionSet(0, 0, -1, 0, 0, 0, 0),Omega_)
       );

       Info << Omega << endl;       




    // PROPERTIES CONCERNING GATHERING STATISTICS

       // Gather/write statistics?
       bool statisticsOn(ABLProperties.lookupOrDefault<bool>("statisticsOn", false));

       // Statistics gathering/writing frequency?
       int statisticsFreq(int(readScalar(ABLProperties.lookup("statisticsFrequency"))));
