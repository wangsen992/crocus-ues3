#include "nesting_utils.H"
#include <fstream>

std::shared_ptr<OGRSpatialReference> ReadProj(const std::string& filename)
{
  std::ifstream proj_file(filename);
  std::string prj;
  std::getline(proj_file, prj);
  proj_file.close();
  std::unique_ptr<OGRSpatialReference> pSRS(new OGRSpatialReference);
  pSRS->SetProjCS(filename.c_str());
  std::cout << prj << std::endl;
  pSRS->SetFromUserInput(prj.c_str());

  return std::move(pSRS);
}

void printSRS(std::shared_ptr<OGRSpatialReference> psrs)
{
  char *pszWKT = nullptr;
  const char* apszOptions[] = { "FORMAT=WKT2_2018", "MULTILINE=YES", nullptr };
  psrs->exportToWkt(&pszWKT, apszOptions);
  printf("%s\n", pszWKT);
  CPLFree(pszWKT);
  
}

void getProjAtts(const netCDF::NcFile& f, WRF_PROJ_PARAMS& params)
{
  f.getAtt("MAP_PROJ").getValues(&params.MAP_PROJ);
  f.getAtt("TRUELAT1").getValues(&params.TRUELAT1);
  f.getAtt("TRUELAT2").getValues(&params.TRUELAT2);
  f.getAtt("MOAD_CEN_LAT").getValues(&params.MOAD_CEN_LAT);
  f.getAtt("STAND_LON").getValues(&params.STAND_LON);
  f.getAtt("POLE_LAT").getValues(&params.POLE_LAT);
  f.getAtt("POLE_LON").getValues(&params.POLE_LON);
}

std::shared_ptr<OGRSpatialReference> getCRS(const WRF_PROJ_PARAMS& params)
{
  std::shared_ptr<OGRSpatialReference> pcrs(new OGRSpatialReference);
  pcrs->SetLCC
  (
    params.TRUELAT1,
    params.TRUELAT2,
    params.MOAD_CEN_LAT,
    params.STAND_LON,
    0,
    0
  );
  return pcrs;

}

std::shared_ptr<OGRSpatialReference> getGCS(std::string gcs_str)
{
  std::shared_ptr<OGRSpatialReference> pgcs(new OGRSpatialReference);
  pgcs->SetWellKnownGeogCS(gcs_str.c_str());
  return pgcs;
}

