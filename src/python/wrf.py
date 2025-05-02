import numpy as np
import xarray as xr

if __name__ == "__main__":
    
    fpath = "../../wrfout_d03_2018-06-25_00_00_00"

    ds = xr.load_dataset(fpath)
    U = ds["U"]
    V = ds["V"]
    PH = ds["PH"]
    PHB = ds["PHB"]
    z = (PH + PHB) / 9.8
