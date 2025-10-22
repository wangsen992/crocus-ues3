import geojson
import argparse, sys
from pyproj import Transformer
from pyproj import CRS
import numpy as np

if __name__ == "__main__":
    parser = \
        argparse.ArgumentParser(prog="get_las_bnbs.py",
                                usage="get argument for keep_xy parameter for las2las ",
                                description=None)

    parser.add_argument("--proj_fname", "-p", default="proj4str.txt")
    parser.add_argument("--subset_geojson")
    parser.add_argument("--las_proj_fname", default="las_proj.txt")

    try:
        args = parser.parse_args(sys.argv[1:])
    except:
        parser.print_help()
        exit()

    # target crs
    with open(args.proj_fname, "r") as f:
        proj_str = f.read()
    target_crs = CRS.from_proj4(proj_str)

    # Las file CRS
    with open(args.las_proj_fname, 'r') as f:
        s = f.read()
    source_crs = CRS.from_wkt(s)

    proj = Transformer.from_crs(source_crs, target_crs)

    # load influence region
    with open(args.subset_geojson,'r') as f:
        gs = geojson.load(f)
    
    coords = np.array(gs["features"][0]["geometry"]["coordinates"]).squeeze()
    x1, y1 = coords.max(axis=0)
    x0, y0 = coords.min(axis=0)

    ((lasx0, lasx1), (lasy0, lasy1)) = proj.transform((x0,x1),(y0,y1),
                                                      direction="INVERSE")

    print(f"{int(lasx0)} {int(lasy0)} {int(lasx1)} {int(lasy1)}")


