"""create boundaries for city4cfd and cfd domain in geojson form"""

from pyproj import Transformer
from pyproj import CRS
import geojson
from geojson import Feature, MultiPolygon, FeatureCollection
import argparse
import sys

def create_bnd(x0, y0, buffer, name):
    mp = MultiPolygon([
        ([(x0-buffer, y0-buffer), (x0+buffer, y0-buffer), 
         (x0+buffer, y0+buffer), (x0-buffer, y0+buffer),])
    ])
    mp['coordinates']= [mp['coordinates']]
    
    feature = Feature(geometry=mp)
    
    feature_collection = FeatureCollection([feature])
    
    feature_collection['name'] = name
    return feature_collection

if __name__ == "__main__":


    ## Input
    parser = \
        argparse.ArgumentParser(prog="create_bnds.py",
                                usage="Create boundaries for city4cfd and "
                                      "cfd domain in geojson form",
                                description=None)

    # proj string 
    parser.add_argument("--proj_fname", "-p", default="proj4str.txt")
    parser.add_argument("--lon0", type=float)
    parser.add_argument("--lat0", type=float)
    parser.add_argument("--building_buffer", type=int, default=200)
    parser.add_argument("--domain_buffer", type=int,
                        default=300)
    parser.add_argument("--target_dir", default="./results")
    try:
        args = parser.parse_args(sys.argv[1:])
    except:
        parser.print_help()
        exit()

    proj_fname = args.proj_fname
    # # center in latlon 
    lon0, lat0 = args.lon0, args.lat0
    # # build and domain size
    c4c_buf = args.building_buffer  # buffer size for city4cfd building
    domain_buf = args.domain_buffer # buffer size for cfd domain

    print(args)


    with open(proj_fname, "r") as f:
        proj_str = f.read()

    # transformation from las-crs to geodetic
    proj_crs = CRS.from_proj4(proj_str)
    proj_st = Transformer.from_crs(proj_crs, proj_crs.geodetic_crs)

    x0, y0 = proj_st.transform(lon0, lat0, direction="INVERSE")
    print(x0, y0)

    c4c_bnd = create_bnd(x0,y0, c4c_buf, "influenceRegion")
    domain_bnd = create_bnd(x0, y0, domain_buf, "domainBnd")

    with open(f"{args.target_dir}/influenceRegion.geojson", "w") as f:
        f.write(geojson.dumps(c4c_bnd))
    with open(f"{args.target_dir}/domainBnd.geojson", "w") as f:
        f.write(geojson.dumps(domain_bnd))
