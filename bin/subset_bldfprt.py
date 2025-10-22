import geojson
import geopandas as gpd
import numpy as np
from pyproj import CRS,Transformer
import argparse, sys

if __name__ == "__main__":

    ## Input
    parser = \
        argparse.ArgumentParser(prog="subset_bldfprt.py",
                                usage="subset the building geojson and "
                                      "transform to target projection",
                                description=None)
    parser.add_argument("--proj_fname", "-p", default="proj4str.txt")
    parser.add_argument("--source_bldfprt")
    parser.add_argument("--target_dir", default="./results")

    try:
        args = parser.parse_args(sys.argv[1:])
    except:
        parser.print_help()
        exit()
    
    with open(args.proj_fname, "r") as f:
        proj_str = f.read()
    
    crs = CRS.from_proj4(proj_str)
    proj_st = Transformer.from_crs(crs, crs.geodetic_crs)
    
    with open(f"{args.target_dir}/influenceRegion.geojson",'r') as f:
        gs = geojson.load(f)
    
    coords = np.array(gs["features"][0]["geometry"]["coordinates"]).squeeze() 
    x1, y1 = coords.max(axis=0)
    x0, y0 = coords.min(axis=0)
    lon0, lat0 = proj_st.transform(x0, y0)
    lon1, lat1 = proj_st.transform(x1, y1)

    
    # This can be optimized by using dask for parallel processing
    # gdf = gpd.read_file("../DATA/Buildings_20250304.csv")
    gdf = gpd.read_file(args.source_bldfprt)
    # Solution from : https://stackoverflow.com/questions/61122875/geopandas-how-to-read-a-csv-and-convert-to-a-geopandas-dataframe-with-polygons
    gdf = gpd.GeoDataFrame(
        geometry=gpd.GeoSeries.from_wkt(gdf['the_geom'], crs='latlon'), data=gdf
    )
    gdf = gdf.set_crs("latlon")
    
    gdf_subset = gdf.cx[lon0:lon1, lat0:lat1]
    
    gdf_subset = gdf_subset.to_crs(crs)
    
    gdf_subset.to_file(f"{args.target_dir}/buildings.geojson", driver="GeoJSON")
