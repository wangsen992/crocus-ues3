import time
from collections import Counter
import io

import laspy
import numpy as np
from itertools import product, combinations
from scipy.spatial import cKDTree
import pyvista as pv

import multiprocessing
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import argparse, sys

def voxel(ar, bottom, vwidth):
    # xb, yb = bottom[:,0], bottom[:,1]

    domain_min = np.min([ar.min(axis=0), bottom.min(axis=0)], axis=0)
    ar_i = ((ar - domain_min)/vwidth).astype("int")
    bottom_i = ((bottom - domain_min)/vwidth).astype("int")
    domain_i_max = np.max([ar_i.max(axis=0), bottom_i.max(axis=0)], axis=0)
    xi_counter = np.unique(ar_i, axis=0, return_counts=True)
    xbi_counter = np.unique(bottom_i, axis=0, return_counts=True)

    vcnt = np.zeros(domain_i_max+1)
    vbcnt = np.zeros(domain_i_max[:2]+1)
    
    # basically counting the number of each unique 3-tuple (xi,yi,zi)
    xi_tuple = xi_counter[0]
    xi_counts = xi_counter[1]
    for i in range(xi_counts.size):
        t = xi_tuple[i]
        vcnt[t[0]-1, t[1]-1, t[2]-1] = xi_counts[i]

    xbi_tuple = xbi_counter[0]
    xbi_counts = xbi_counter[1]
    for i in range(xbi_counts.size):
        t = xbi_tuple[i]
        if t[0] >= 0 and t[1] >= 0:
            vbcnt[t[0]-1, t[1]-1] = xbi_counts[i]

    x,y,z = ar[:,0], ar[:,1], ar[:,2]
    vx = np.arange(x.min(), x.max()+vwidth, vwidth)
    vy = np.arange(y.min(), y.max()+vwidth, vwidth)
    vz = np.arange(z.min(), z.max()+vwidth, vwidth)
    return vx, vy, vz, vcnt, vbcnt

def to_OF_list(name, list_data):
    return f"{name} nonuniform {len(list_data)}" + "\n" + "(\n\t" + "\n\t".join(list_data) + "\n);"

def proc_subset(fname_veg):
    veg_las = laspy.read(fname_veg)

    if len(veg_las) < 10:
        print(f"No veg points {fname_veg}")
        return [], []
    subset_prefix = fname_veg.stem.split("_veg")[0]
    fname_gnd = gnd_subset_dir/f"{subset_prefix}_ground.las"
    ground_las = laspy.read(fname_gnd)

    print(f"Processing {fname_veg}")

    # extract valid points from veg and ground
    veg_points = veg_las.xyz
    veg_scan_angle = 0.006 * veg_las.scan_angle
    terrain_points = ground_las.xyz
    terrain_scan_angle = 0.006 * ground_las.scan_angle

    # voxelise 
    spacing = 0.5
    x, y, z, vcnt, vbcnt = voxel(veg_points, terrain_points, spacing)

    # Pyvista stuff
    cond = np.abs(terrain_scan_angle) < 2
    terrain_point_cloud = pv.PolyData(terrain_points[cond, :])
    terrain_point_cloud.point_data["scan_angle"] = terrain_scan_angle[cond]

    veg_mesh = pv.PolyData(veg_points)
    veg_mesh.point_data["height"] = veg_points[:, 2]
    veg_mesh.point_data["scan_angle"] = veg_scan_angle

    grid = pv.ImageData()
    grid.dimensions = np.array(vcnt.shape)+1
    grid.origin = veg_points.min(axis=0)
    grid.spacing = np.ones(3) * spacing

    # Compute LAD with Beer-Lambert theory based method
    vcnt_cb = np.concat([vbcnt[...,np.newaxis], vcnt], axis=2)
    vcnt_cb_sum = vcnt_cb.cumsum(axis=2)
    vcnt_cb_sum[vcnt_cb_sum == 0] = 1

    del vcnt_cb, vbcnt, vcnt
    k = 0.5
    vcnt_cb_sum_tmp = np.log(vcnt_cb_sum[:,:,1:]/ vcnt_cb_sum[:,:,:-1])

    lad = vcnt_cb_sum_tmp/ (k * spacing)

    # Output
    points_list = []
    lad_list = []

    for i in range(x.size-1):
        for j in range(y.size-1):
            for k in range(z.size-1):
                if lad[i,j,k] > 0:
                    points_list.append(f"({x[i]:6.2f} {y[j]:6.2f} {z[k]:6.2f})")
                    lad_list.append(f"{lad[i,j,k]:4.2f}")

    return points_list, lad_list

def to_foam(points_list, lad_list):
    header = """
    FoamFile
    {
        verison 2.0;
        format  ascii;
        class   dictionary;
        object  lad;
    }
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    """

    points_entry = to_OF_list("points", points_list)
    lad_entry = to_OF_list("lad", lad_list)
    dimension_entry = "dimension  [0 -1 0 0 0];\n"
    print("Write data as OpenFOAM dicts directly...")
    with open("constant/urban/point_data/lad",'w', buffering=8192) as f:
        f.write(header)
        f.write(dimension_entry)
        print("header written")
        f.write(points_entry)
        print("points written")
        f.write(lad_entry)
        print("lad written")

if __name__ == "__main__":
    ## Input
    parser = \
        argparse.ArgumentParser(prog="canopy_las2of.py",
                                usage="convert las to openfoam lad list ",
                                description=None)
    parser.add_argument("--num_workers", default="4")
    try:
        args = parser.parse_args(sys.argv[1:])
    except:
        parser.print_help()
        exit()
    num_workers = int(args.num_workers)
    # read ground and vegetation points
    veg_subset_dir = Path("ppcfd_results/vegetation_las")
    gnd_subset_dir = Path("ppcfd_results/ground_las")

    points_list, lad_list = [], []
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        for pl, ll in executor.map(proc_subset, list(veg_subset_dir.iterdir())):
            points_list += pl
            lad_list += ll

    # to_foam(points_list, lad_list)
