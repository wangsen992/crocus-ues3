#include "nesting_utils.H"
#include "emptyPolyPatch.H"
#include "cellModeller.H"


void readWrfCaseInfo(struct WrfCaseInfo * wrfInfo, netCDF::NcFile& dataFile)
{
    
      wrfInfo->DateStrLen = dataFile.getDim("DateStrLen").getSize();
      wrfInfo->N_Time = dataFile.getDim("Time").getSize();
      wrfInfo->N_west_east = dataFile.getDim("west_east").getSize();
      wrfInfo->N_west_east_stag = dataFile.getDim("west_east_stag").getSize();
      wrfInfo->N_south_north = dataFile.getDim("south_north").getSize();
      wrfInfo->N_south_north_stag = dataFile.getDim("south_north_stag").getSize();
      wrfInfo->N_bottom_top = dataFile.getDim("bottom_top").getSize();
      wrfInfo->N_bottom_top_stag = dataFile.getDim("bottom_top_stag").getSize();

      dataFile.getAtt("MAP_PROJ").getValues(&(wrfInfo->MAP_PROJ));
      dataFile.getAtt("CEN_LAT").getValues(&(wrfInfo->CEN_LAT));
      dataFile.getAtt("CEN_LON").getValues(&(wrfInfo->CEN_LON));
      dataFile.getAtt("TRUELAT1").getValues(&(wrfInfo->TRUELAT1));
      dataFile.getAtt("TRUELAT2").getValues(&(wrfInfo->TRUELAT2));
      dataFile.getAtt("MOAD_CEN_LAT").getValues(&(wrfInfo->MOAD_CEN_LAT));
      dataFile.getAtt("STAND_LON").getValues(&(wrfInfo->STAND_LON));
      dataFile.getAtt("POLE_LAT").getValues(&(wrfInfo->POLE_LAT));
      dataFile.getAtt("POLE_LON").getValues(&(wrfInfo->POLE_LON));
      dataFile.getAtt("DX").getValues(&(wrfInfo->DX));
      dataFile.getAtt("DY").getValues(&(wrfInfo->DY));

      wrfInfo->Nx = wrfInfo->N_west_east_stag-2;
      wrfInfo->Ny = wrfInfo->N_south_north_stag-2;
      wrfInfo->Nz = wrfInfo->N_bottom_top_stag-2;
      wrfInfo->Ncellx = wrfInfo->Nx-1;
      wrfInfo->Ncelly = wrfInfo->Ny-1;
      wrfInfo->Ncellz = wrfInfo->Nz-1;

      wrfInfo->timestamps.resize(wrfInfo->N_Time);
      char* tmp_time = new char[wrfInfo->N_Time*wrfInfo->DateStrLen];
      netCDF::NcVar Times(dataFile.getVar("Times"));
      Times.getVar(tmp_time);

      for (int i=0; i < wrfInfo->N_Time; i++)
      {
        std::string& si = wrfInfo->timestamps[i];
        for (int j=0; j < wrfInfo->DateStrLen; j++)
        {
          si.push_back( *(tmp_time + i * wrfInfo->DateStrLen + j)); 
        }
      }

}

void testFunc(netCDF::NcFile& dataFile)
{
}

using namespace Foam;
autoPtr<polyMesh> meshFromNc(netCDF::NcFile& dataFile, Time& runTime, bool write=false)
{
    // Start parsing the coordinate systems of the result
      // Collect dimension sizes
      int N_Time(dataFile.getDim("Time").getSize());
      int N_west_east(dataFile.getDim("west_east").getSize());
      int N_west_east_stag(dataFile.getDim("west_east_stag").getSize());
      int N_south_north(dataFile.getDim("south_north").getSize());
      int N_south_north_stag(dataFile.getDim("south_north_stag").getSize());
      int N_bottom_top(dataFile.getDim("bottom_top").getSize());
      int N_bottom_top_stag(dataFile.getDim("bottom_top_stag").getSize());

      int Nx(N_west_east_stag-2);
      int Ny(N_south_north_stag-2);
      int Nz(N_bottom_top_stag-2);
      int Ncellx(Nx-1);
      int Ncelly(Ny-1);
      int Ncellz(Nz-1);

      std::shared_ptr<OGRSpatialReference> pgcs(new OGRSpatialReference);
      pgcs->SetWellKnownGeogCS("WGS84");
      WRF_PROJ_PARAMS params;
      getProjAtts(dataFile, params);
      std::shared_ptr<OGRSpatialReference> pcrs=getCRS(params);
      std::shared_ptr<OGRCoordinateTransformation> poCT(OGRCreateCoordinateTransformation(pgcs.get(), pcrs.get()));
      
      struct WrfCaseInfo wrfInfo;
      readWrfCaseInfo(&wrfInfo, dataFile);

      autoPtr<polyMesh> tmesh;
      
      // Loading coordinates from lat/lon to xy
      float tmp_full[N_Time][N_south_north_stag][N_west_east_stag];
      // Load lat and lon (note only first time step is needed)
      float X[N_south_north][N_west_east];
      float Y[N_south_north][N_west_east];

      {
        float xlat[N_south_north][N_west_east];
        netCDF::NcVar XLAT(dataFile.getVar("XLAT"));
        XLAT.getVar(tmp_full);
        for(int isn=0; isn < N_south_north; isn++)
        {
            for(int iwe=0; iwe < N_west_east; iwe++)
            {
                xlat[isn][iwe] = tmp_full[0][isn][iwe];
            }
        }

        float xlon[N_south_north][N_west_east];
        netCDF::NcVar XLON(dataFile.getVar("XLONG"));
        XLON.getVar(tmp_full);
        for(int isn=0; isn < N_south_north; isn++)
        {
            for(int iwe=0; iwe < N_west_east; iwe++)
            {
                xlon[isn][iwe] = tmp_full[0][isn][iwe];
            }
        }
        
        double x,y;
        for(int isn=0; isn < N_south_north; isn++)
        {
            for(int iwe=0; iwe < N_west_east; iwe++)
            {
                x = xlat[isn][iwe];
                y = xlon[isn][iwe];
                poCT->Transform(1, &x, &y);
                X[isn][iwe] = x;
                Y[isn][iwe] = y;

                // std::cout << "X, Y" << ": " << X[isn][iwe] << ", " << Y[isn][iwe] << std::endl;
            }
        }
      }

      // Load lat and lon (U_stag, note only first time step is needed)
      float X_U[N_south_north][N_west_east_stag];
      float Y_U[N_south_north][N_west_east_stag];
      {
        float xlat_u[N_south_north][N_west_east_stag];
        // float tmp_full_U_stag[N_Time][N_south_north][N_west_east_stag];
        netCDF::NcVar XLAT_U(dataFile.getVar("XLAT_U"));
        XLAT_U.getVar(tmp_full);
        for(int isn=0; isn < N_south_north; isn++)
        {
            for(int iwe=0; iwe < N_west_east_stag; iwe++)
            {
                xlat_u[isn][iwe] = tmp_full[0][isn][iwe];
            }
        }

        float xlon_u[N_south_north][N_west_east_stag];
        netCDF::NcVar XLON_U(dataFile.getVar("XLONG_U"));
        XLON_U.getVar(tmp_full);
        for(int isn=0; isn < N_south_north; isn++)
        {
            for(int iwe=0; iwe < N_west_east_stag; iwe++)
            {
                xlon_u[isn][iwe] = tmp_full[0][isn][iwe];
            }
        }

        double x,y;
        for(int isn=0; isn < N_south_north; isn++)
        {
            for(int iwe=0; iwe < N_west_east_stag; iwe++)
            {
                x = xlat_u[isn][iwe];
                y = xlon_u[isn][iwe];
                poCT->Transform(1, &x, &y);
                X_U[isn][iwe] = x;
                Y_U[isn][iwe] = y;
                // std::cout << "X, Y" << ": " << X[isn][iwe] << ", " << Y[isn][iwe] << std::endl;
            }
        }
      }
      
      // Load lat and lon (V_stag, note only first time step is needed)
      // float tmp_full_V_stag[N_Time][N_south_north_stag][N_west_east];
      float X_V[N_south_north_stag][N_west_east];
      float Y_V[N_south_north_stag][N_west_east];
      {
        float xlat_v[N_south_north_stag][N_west_east];
        netCDF::NcVar XLAT_V(dataFile.getVar("XLAT_V"));
        XLAT_V.getVar(tmp_full);
        for(int isn=0; isn < N_south_north_stag; isn++)
        {
            for(int iwe=0; iwe < N_west_east; iwe++)
            {
                xlat_v[isn][iwe] = tmp_full[0][isn][iwe];
            }
        }

        float xlon_v[N_south_north_stag][N_west_east];
        netCDF::NcVar XLON_V(dataFile.getVar("XLONG_V"));
        XLON_V.getVar(tmp_full);
        for(int isn=0; isn < N_south_north_stag; isn++)
        {
            for(int iwe=0; iwe < N_west_east; iwe++)
            {
                xlon_v[isn][iwe] = tmp_full[0][isn][iwe];
            }
        }

        double x,y;
        for(int isn=0; isn < N_south_north_stag; isn++)
        {
            for(int iwe=0; iwe < N_west_east; iwe++)
            {
                x = xlat_v[isn][iwe];
                y = xlon_v[isn][iwe];
                poCT->Transform(1, &x, &y);
                X_V[isn][iwe] = x;
                Y_V[isn][iwe] = y;
            }
        }
      }
      
      // Load z (need to assign real value)
      int cc = 0;
      float * Z_stag = new float[N_bottom_top_stag * N_south_north * N_west_east];
      {
        float* phb = new float[N_Time*N_bottom_top_stag*N_south_north*N_west_east];
        float* ph = new float[N_Time*N_bottom_top_stag*N_south_north*N_west_east];
        netCDF::NcVar phb_nc(dataFile.getVar("PHB"));
        netCDF::NcVar ph_nc(dataFile.getVar("PH"));
        phb_nc.getVar(phb);
        ph_nc.getVar(ph);
        double Re(6371000.0);
        double Phi;

        for (int k=0; k < N_bottom_top_stag; k++)
        {
          for (int j=0; j < N_south_north; j++)
          {
            for (int i=0; i < N_west_east; i++)
            {
              cc = k * N_south_north * N_west_east + j * N_west_east + i;
              // Z_stag[cc] = (phb[cc] + phb[cc]) / 9.81;
              Phi = ph[cc] + phb[cc];
              Z_stag[cc] = (Phi * Re) / (9.81 * Re - Phi);
            }
          }
        }
      }
     
      // Construct OpenFOAM mesh from the given xy points
      // Note: the outer most columes and rows are ignored, and 
      //       the mesh is slightly distored due to two stagged grid are
      //       combined directly
       
        // Construct points using stagged grid points
        int N_points(Nx*Ny*Nz);
        pointField points(N_points);
        int pc=0;
        for(int i=1; i < N_west_east_stag-1; i++)
        {
          for(int j=1; j < N_south_north_stag-1; j++)
          {
            for(int k=1; k < N_bottom_top_stag-1; k++)
            {
                pc = (i-1) + Nx*(j-1) + Nx*Ny*(k-1);
                cc = k * N_south_north * N_west_east + j * N_west_east + i;
                points[pc] = Foam::vector(X_U[j][i], Y_V[j][i], Z_stag[cc]);
            }
          }
        }
        
        // Construct Cells
        List<FixedList<label, 8>> cells(Ncellx*Ncelly*Ncellz);
        for(int i=0; i < Ncellx; i++)
        {
          for(int j=0; j < Ncelly; j++)
          {
            for(int k=0; k < Ncellz; k++)
            {
                cc= i + Ncellx*j + Ncellx*Ncelly*k;
                cells[cc][0]= i + Nx*j + Nx*Ny*k;
                cells[cc][1]= (i+1) + Nx*j + Nx*Ny*k;
                cells[cc][2]= (i+1) + Nx*(j+1) + Nx*Ny*k;
                cells[cc][3]= i + Nx*(j+1) + Nx*Ny*k;
                cells[cc][4]= i + Nx*j + Nx*Ny*(k+1);
                cells[cc][5]= (i+1) + Nx*j + Nx*Ny*(k+1);
                cells[cc][6]= (i+1) + Nx*(j+1) + Nx*Ny*(k+1);
                cells[cc][7]= i + Nx*(j+1) + Nx*Ny*(k+1);
            }
          }
        }
        const cellModel& hex = *(cellModeller::lookup("hex"));
        cellShapeList cellShapes(cells.size());
        for (int celli = 0; celli < cells.size(); celli++)
        {
            cellShapes[celli] = cellShape(hex, labelList(cells[celli]), true);
        }

        // Now set the patches (borrwo from createBlock.C)
        FixedList<List<FixedList<label, 4>>, 6> patches;
        label patchi = 0;
        label facei = 0;
        // x-direction

        // x-min
        patches[patchi].setSize(Ncelly*Ncellz);
        for (label k=0; k<Ncellz; k++)
        {
            for (label j=0; j<Ncelly; j++)
            {
                patches[patchi][facei][0] = 0 + Nx*j + Nx*Ny*k;
                patches[patchi][facei][1] = 0 + Nx*j + Nx*Ny*(k+1);
                patches[patchi][facei][2] = 0 + Nx*(j+1) + Nx*Ny*(k+1);
                patches[patchi][facei][3] = 0 + Nx*(j+1) + Nx*Ny*k;


                facei++;
            }
        }

        // x-max
        patchi++;
        facei = 0;

        patches[patchi].setSize(Ncelly*Ncellz);

        for (label k=0; k<Ncellz; k++)
        {
            for (label j=0; j<Ncelly; j++)
            {
                patches[patchi][facei][0] = Ncellx + Nx*j + Nx*Ny*k;
                patches[patchi][facei][1] = Ncellx + Nx*(j+1) + Nx*Ny*k;
                patches[patchi][facei][2] = Ncellx + Nx*(j+1) + Nx*Ny*(k+1);
                patches[patchi][facei][3] = Ncellx + Nx*j + Nx*Ny*(k+1);
                facei++;
            }
        }

        // y-direction

        // y-min
        patchi++;
        facei = 0;

        patches[patchi].setSize(Ncellx*Ncellz);
        for (label i=0; i<Ncellx; i++)
        {
            for (label k=0; k<Ncellz; k++)
            {
                patches[patchi][facei][0] = i + 0 + Nx*Ny*k;
                patches[patchi][facei][1] = i+1 + 0 + Nx*Ny*k;
                patches[patchi][facei][2] = i+1 + 0 + Nx*Ny*(k+1);
                patches[patchi][facei][3] = i + 0 + Nx*Ny*(k+1);
                facei++;
            }
        }

        // y-max
        patchi++;
        facei = 0;

        patches[patchi].setSize(Ncellx*Ncellz);

        for (label i=0; i<Ncellx; i++)
        {
            for (label k=0; k<Ncellz; k++)
            {
                patches[patchi][facei][0] = i   + Nx*Ncelly + Nx*Ny*k;
                patches[patchi][facei][1] = i + Nx*Ncelly + Nx*Ny*(k+1);
                patches[patchi][facei][2] = i+1 + Nx*Ncelly + Nx*Ny*(k+1);
                patches[patchi][facei][3] = i+1   + Nx*Ncelly + Nx*Ny*k;
                facei++;
            }
        }

        // z-direction

        // z-min
        patchi++;
        facei = 0;

        patches[patchi].setSize(Ncellx*Ncelly);

        for (label i=0; i<Ncellx; i++)
        {
            for (label j=0; j<Ncelly; j++)
            {
                patches[patchi][facei][0] = i + Nx*j + 0;
                patches[patchi][facei][1] = i + Nx*(j+1) + 0;
                patches[patchi][facei][2] = (i+1) + Nx*(j+1) + 0;
                patches[patchi][facei][3] = (i+1) + Nx*j + 0;
                facei++;
            }
        }

        // z-max
        patchi++;
        facei = 0;

        patches[patchi].setSize(Ncellx*Ncelly);

        for (label i=0; i<Ncellx; i++)
        {
            for (label j=0; j<Ncelly; j++)
            {
                patches[patchi][facei][0] = i   + Nx*j     + Nx*Ny*Ncellz;
                patches[patchi][facei][1] = (i+1)   + Nx*j + Nx*Ny*Ncellz;
                patches[patchi][facei][2] = (i+1) + Nx*(j+1) + Nx*Ny*Ncellz;
                patches[patchi][facei][3] = i + Nx*(j+1)     + Nx*Ny*Ncellz;

                facei++;
            }
        }

        faceListList patchLists(6);
        forAll(patches, patchi)
        {
          patchLists[patchi].setSize(patches[patchi].size());
          forAll(patches[patchi], facei)
          {
            patchLists[patchi][facei] = face(labelList(patches[patchi][facei]));
          }
        }

        wordList patchNames(6);
        patchNames[0] = "west";
        patchNames[1] = "east";
        patchNames[2] = "south";
        patchNames[3] = "north";
        patchNames[4] = "bottom";
        patchNames[5] = "top";

        PtrList<dictionary> patchDicts(6);
        for (int i = 0; i < 6; i++)
        {
            patchDicts.set(i, new dictionary());
            patchDicts[i].add("type", "zeroGradient");
        }
        Info << "Constructing mesh." << endl;

        // Construct polyMesh
        tmesh.set
        (
          new polyMesh
          (
            IOobject
            (
              polyMesh::defaultRegion,
              runTime.constant(),
              runTime
            ),
            std::move(points),
            cellShapes,
            patchLists,
            patchNames,
            patchDicts,
            "defaultFaces",
            emptyPolyPatch::typeName
          )
        );
    if(write)
    {
      tmesh->write();
    }

    return tmesh;
}

Foam::autoPtr<Foam::fvMesh> fvmeshFromNc(netCDF::NcFile& ncfile, Foam::Time& runTime)
{
  Foam::autoPtr<Foam::polyMesh> polyMesh
  (
    meshFromNc(ncfile, runTime, true)
  );
  Foam::autoPtr<Foam::fvMesh> fvmesh
  (
    new Foam::fvMesh
    (
      IOobject
      (
        fvMesh::defaultRegion,
        runTime.timeName(),
        runTime
      )
    )
  );
  Info << "fvMesh completed." << endl;
  return fvmesh;
}
    

