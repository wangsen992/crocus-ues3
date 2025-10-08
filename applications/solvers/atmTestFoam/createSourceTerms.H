    // Create the source term, which lumps together the driving pressure gradient and mesoscale source terms
    // for advection of momentum and potential temperature.

    // Create the momentum sourc term.
    Info << "Creating the momentum source term, SourceU..." << endl;
    volVectorField SourceU
    (
        IOobject
        (
            "SourceU",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("SourceU",dimVelocity/dimTime,vector::zero)
    );
    volVectorField Source_dUdt_advec
    (
        IOobject
        (
            "Source_dUdt_advec",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("Source_dUdt_advec",dimVelocity/dimTime,vector::zero)
    );
    volVectorField Source_Ug
    (
        IOobject
        (
            "Source_Ug",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("Source_Ug",dimVelocity,vector::zero)
    );
    volVectorField Source_U_LS
    (
        IOobject
        (
            "Source_U_LS",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("Source_U_LS",dimVelocity,vector::zero)
    );

    Info << "Creating the potential temperature source term, SourceT..." << endl;
    volScalarField SourceT
    (
        IOobject
        (
            "SourceT",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("SourceT",dimTemperature/dimTime,0.0)
    );
    volScalarField Source_thl_LS
    (
        IOobject
        (
            "Source_thl_LS",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("Source_thl_LS",dimTemperature,0.0)
    );
    volScalarField Source_dthldt_advec
    (
        IOobject
        (
            "Source_dthldt_advec",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("Source_dthldt_advec",dimTemperature/dimTime,0.0)
    );

    volScalarField Source_qt_LS
    (
        IOobject
        (
            "Source_qt_LS",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("Source_qt_LS",dimless,0.0)
    );
    volScalarField Source_dqtdt_advec
    (
        IOobject
        (
            "Source_dqtdt_advec",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("Source_dqtdt_advec",dimless/dimTime,0.0)
    );

    List<scalar> sourceUXColumn(hLevelsTotal,0.0);
    List<scalar> sourceUYColumn(hLevelsTotal,0.0);
    List<scalar> sourceUZColumn(hLevelsTotal,0.0);
    List<scalar> sourceTColumn(hLevelsTotal,0.0);
    List<scalar> source_dudt_advecColumn(hLevelsTotal,0.0);
    List<scalar> source_dvdt_advecColumn(hLevelsTotal,0.0);
    List<scalar> source_dthldt_advecColumn(hLevelsTotal,0.0);
    List<scalar> source_dqtdt_advecColumn(hLevelsTotal,0.0);
    List<scalar> source_uColumn(hLevelsTotal,0.0);
    List<scalar> source_vColumn(hLevelsTotal,0.0);
    List<scalar> source_wColumn(hLevelsTotal,0.0);
    List<scalar> source_thlColumn(hLevelsTotal,0.0);
    List<scalar> source_qtColumn(hLevelsTotal,0.0);
    List<scalar> source_ugColumn(hLevelsTotal,0.0);
    List<scalar> source_vgColumn(hLevelsTotal,0.0);



    // Create the momentum source term time history directory
    if (Pstream::master() && !isDir(postProcessingDir/"SourceHistory"))
    {
        mkDir(postProcessingDir/"SourceHistory");
    }
    if (Pstream::master() && !isDir(postProcessingDir/"SourceHistory"/runTime.timeName()))
    {
        mkDir(postProcessingDir/"SourceHistory"/runTime.timeName());
    }



    // Create the momentum source term time history files.
    autoPtr<OFstream> sourceUXHistoryFile;
    if (Pstream::master())
    {
        sourceUXHistoryFile.reset(new OFstream(postProcessingDir/"SourceHistory"/runTime.timeName()/"SourceUXHistory"));

        if (nSourceMomentumHeights > 1)
        {
            sourceUXHistoryFile() << "Heights (m) ";
            forAll(hLevelsValues,i)
            {
                sourceUXHistoryFile() << hLevelsValues[i] << " ";
            }
            sourceUXHistoryFile() << endl;
        }

        sourceUXHistoryFile() << "Time (s)" << " " << "dt (s)" << " " << "source term (m/s^2)" << endl;
    }

    autoPtr<OFstream> sourceUYHistoryFile;
    if (Pstream::master())
    {
        sourceUYHistoryFile.reset(new OFstream(postProcessingDir/"SourceHistory"/runTime.timeName()/"SourceUYHistory"));

        if (nSourceMomentumHeights > 1)
        {
            sourceUYHistoryFile() << "Heights (m) ";
            forAll(hLevelsValues,i)
            {
                sourceUYHistoryFile() << hLevelsValues[i] << " ";
            }
            sourceUYHistoryFile() << endl;
        }

        sourceUYHistoryFile() << "Time (s)" << " " << "dt (s)" << " " << "source term (m/s^2)" << endl;
    }

    autoPtr<OFstream> sourceUZHistoryFile;
    if (Pstream::master())
    {
        sourceUZHistoryFile.reset(new OFstream(postProcessingDir/"SourceHistory"/runTime.timeName()/"SourceUZHistory"));

        if (nSourceMomentumHeights > 1)
        {
            sourceUZHistoryFile() << "Heights (m) ";
            forAll(hLevelsValues,i)
            {
                sourceUZHistoryFile() << hLevelsValues[i] << " ";
            }
            sourceUZHistoryFile() << endl;
        }

        sourceUZHistoryFile() << "Time (s)" << " " << "dt (s)" << " " << "source term (m/s^2)" << endl;
    }



    // Create the potential temperature source term time history files.
    autoPtr<OFstream> sourceTHistoryFile;
    if (Pstream::master())
    {
        sourceTHistoryFile.reset(new OFstream(postProcessingDir/"SourceHistory"/runTime.timeName()/"SourceTHistory"));

        if (nSourceMomentumHeights > 1)
        {
            sourceTHistoryFile() << "Heights (m) ";
            forAll(hLevelsValues,i)
            {
                sourceTHistoryFile() << hLevelsValues[i] << " ";
            }
            sourceTHistoryFile() << endl;
        }

        sourceTHistoryFile() << "Time (s)" << " " << "dt (s)" << " " << "source term (K/s)" << endl;
    }
