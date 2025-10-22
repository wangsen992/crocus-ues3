"""transform raw las to target proj and separate building and ground"""

import laspy
import numpy as np
from pyproj import Transformer
from pyproj import CRS
import argparse, sys
import geojson
import multiprocessing
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

def proc_las(fname):
    print(fname)

    proj_fname = args.proj_fname
    with open(proj_fname, "r") as f:
        proj_str = f.read()
    
    las = laspy.read(fname)

    # Create transformations
    crs = las.header.parse_crs()

    # transformation from las-crs to tartget projection
    proj = Transformer.from_crs(crs, CRS.from_proj4(proj_str))

    # load influence region
    with open(args.subset_geojson,'r') as f:
        gs = geojson.load(f)

    coords = np.array(gs["features"][0]["geometry"]["coordinates"]).squeeze()
    x1, y1 = coords.max(axis=0)
    x0, y0 = coords.min(axis=0)

    # Subset raw las into target area
    (lasx0, lasx1), (lasy0, lasy1) = proj.transform((x0,x1),(y0,y1), direction='INVERSE')
    
    # Subset raw las into target area
    subset_las = laspy.LasData(las.header)
    cond = (las.x > lasx0) & (las.x < lasx1) & (las.y > lasy0) & (las.y < lasy1)
    subset_las.points = las.points[cond]

    del las


    # Transformation to target crs
    x,y,z = proj.transform(subset_las.x, subset_las.y, subset_las.z)
    z = z * 0.3048
    new_file = laspy.create(point_format=subset_las.point_format, file_version=subset_las.header.version)

    new_file.x = x
    new_file.y = y
    new_file.z = z
    new_file.classification = subset_las.classification
    new_file.scan_angle = subset_las.scan_angle

    new_file_ground = laspy.create(point_format=subset_las.point_format, file_version=subset_las.header.version)
    new_file_building = laspy.create(point_format=subset_las.point_format, file_version=subset_las.header.version)
    new_file_vegetation = laspy.create(point_format=subset_las.point_format, file_version=subset_las.header.version)
    new_file_water = laspy.create(point_format=subset_las.point_format, file_version=subset_las.header.version)
    new_file_ground.points = new_file.points[list(map(lambda x: x in [2,11], new_file.classification))]
    new_file_building.points = new_file.points[new_file.classification==6]
    new_file_vegetation.points = new_file.points[list(map(lambda x: x in [3,4,5], new_file.classification))]
    new_file_water.points = new_file.points[new_file.classification==9]
    del new_file

    new_file_ground.write(ground_dir/f"{fname.stem}_ground.las")
    del new_file_ground
    new_file_building.write(building_dir/f"{fname.stem}_building.las")
    del new_file_building
    new_file_vegetation.write(vegetation_dir/f"{fname.stem}_vegetation.las")
    del new_file_vegetation
    new_file_water.write(water_dir/f"{fname.stem}_water.las")

    return 0

if __name__ == "__main__":
    ## Input
    parser = \
        argparse.ArgumentParser(prog="prep_las.py",
                                usage="transform raw las into target proj "
                                      "and separate building and ground",
                                description=None)
    parser.add_argument("--proj_fname", "-p", default="proj4str.txt")
    parser.add_argument("--source_las")
    parser.add_argument("--subset_geojson")
    parser.add_argument("--target_dir", default="./results")
    parser.add_argument("--num_workers", default="4")
    try:
        args = parser.parse_args(sys.argv[1:])
    except:
        parser.print_help()
        exit()

    subset_las_dir = Path(args.source_las)

    target_dir = Path(args.target_dir)
    ground_dir = target_dir/"ground_las"
    building_dir = target_dir/"building_las"
    vegetation_dir = target_dir/"vegetation_las"
    water_dir = target_dir/"water_las"

    for cat_dir in [ground_dir, building_dir, vegetation_dir, water_dir]:
        cat_dir.mkdir()

    # print(f"subset_las_dir: {subset_las_dir}")
    # for las_name in subset_las_dir.iterdir():
    #     print(las_name)
    #     proc_las(las_name)

    num_workers = int(args.num_workers)
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        results = executor.map(proc_las,
                               list(subset_las_dir.iterdir()))
