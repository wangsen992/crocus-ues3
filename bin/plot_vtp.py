import sys
import numpy as np
import pandas as pd
from pathlib import Path

import pyvista as pv

import time
from tqdm import tqdm

import multiprocessing

pv.set_jupyter_backend('static')

start_time = pd.to_datetime("2024-07-01 00:00:00")

def surf_dir_check(surf_dir):
    """Check the surf_dir structure is as intended

    The structure should be surface (level 0), time step name (level 1), different surface 
    file names in level 2. 
    """

    return True
    
def surface_file_gen(surf_dir, fname, max_files=None):

    assert surf_dir_check(surf_dir), \
           "input surf_dir must conform to openfoam postProcessing surfaces structure"

    d_len = len(list(surf_dir.iterdir()))
    if max_files is None:
        max_files = d_len
        
    count = 0
    for d in sorted(surf_dir.iterdir(), key=lambda f: float(f.name)):
        if count >= max_files:
            break
        if (d/fname).exists():
            yield float(d.name), pv.read(d/fname)
        else:
            continue
        count += 1

        
def animate(fname, vname, movie_name=None, max_files=None, mesh_kwargs={}, tqdm_pos=0):
    if max_files is None:
        max_files = d_len
    plotter = pv.Plotter()
    plotter.open_movie(movie_name, framerate=30)
    for t, f_z2 in tqdm(surface_file_gen(surface_files_dir, fname, max_files), 
                        total=max_files,
                        desc=movie_name,
                        position=tqdm_pos
                        ):

        time_str = (start_time + pd.to_timedelta(t, unit='s'))\
                            .isoformat(sep=" ", timespec='seconds')

        # Update Z and write a frame for each updated position

        try:
            plotter.add_mesh(f_z2, scalars=vname, **mesh_kwargs)
            plotter.view_xy()
            plotter.add_text(f"time: {time_str} ", 
                             position="upper_edge", 
                             font_size=8, 
                             color='black')
            # Write a frame. This triggers a render.
            plotter.write_frame()
            plotter.clear()
        
        except:
            continue
    # Closes and finalizes movie
    plotter.close()

if __name__ == "__main__":
    print("Offscreen rendering starts....")
    pv.OFF_SCREEN = True

    # Parameters
    surface_files_dir = Path("postProcessing/surfaces/")
    d_len = len(list(surface_files_dir.iterdir()))
    print(d_len)

    # Options
    scalar_bar_args = {'position_x': 0.9, 'position_y': 0.3, 'vertical': True}

    T2_args = {
        "fname": "two_meter_terrain.vtp", 
        "vname": "T", 
        "movie_name": "two_meter_terrain_T.mp4", 
        "mesh_kwargs": {
                "cmap" : "coolwarm",
                "scalar_bar_args": scalar_bar_args
                },
        "tqdm_pos": 0
    }
    U2_args = {
        "fname": "two_meter_terrain.vtp", 
        "vname": "U", 
        "movie_name": "two_meter_terrain_U.mp4", 
        "mesh_kwargs": {
                "cmap" : "jet",
                "clim" : (0, 4.5),
                "scalar_bar_args": scalar_bar_args
                },
        "tqdm_pos": 1
    }
    T10_args = {
        "fname": "ten_meter_terrain.vtp", 
        "vname": "T", 
        "movie_name": "ten_meter_terrain_T.mp4", 
        "mesh_kwargs": {
                "cmap" : "coolwarm",
                "clim" : (288, 295),
                "scalar_bar_args": scalar_bar_args
                },
        "tqdm_pos": 0
    }
    U10_args = {
        "fname": "ten_meter_terrain.vtp", 
        "vname": "U", 
        "movie_name": "ten_meter_terrain_U.mp4", 
        "mesh_kwargs": {
                "cmap" : "jet",
                "clim" : (0, 4.5),
                "scalar_bar_args": scalar_bar_args
                },
        "tqdm_pos": 1
    }
    count = 0
    p_list = []
    max_files = 10000
    for args in [T10_args, U10_args]:
        args['tqdm_pos'] = count
        args['max_files'] = max_files
        p = multiprocessing.Process(target=animate, 
                                    kwargs=args)
        p.start()
        p_list.append(p)
        print(count)
        count += 1
