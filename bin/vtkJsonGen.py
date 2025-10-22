#!/bin/python3

import json
import os
import pathlib
import itertools

def combine_vtk(*vtk_files, out_prefix="output"):
    """Combine multiple vtk files with one field to one vtk file with multiple
    fields

    This is done by finding the line prefixed with FIELD and copy everything
    below to the output file, which is first copied from the first entry. 

    This function returns a filepath as handle for the vtk series file. 
    """

    if (len(vtk_files) < 2):
        print("no need for combining")
        return ""
    
    out_fname = out_prefix+".vtp"
    out_fpath = vtk_files[0].absolute().parent / out_fname

    with open(out_fpath,"wb") as out_f:
        with open(vtk_files[0], "rb") as zero_f:
            lines = zero_f.readlines()
            joined_lines = b''.join(lines)
            out_f.write(joined_lines)

        for i, in_fname in enumerate(vtk_files[1:]):
            if out_prefix in in_fname.as_posix():
                continue
            print(f"Reading {in_fname.name}")
            with open(in_fname, "rb") as in_f:
                w_switch = False
                while True:
                    in_str = in_f.readline()
                    if len(in_str) == 0:
                        break
                    if in_str.startswith(b'FIELD'):
                        w_switch = True
                    if w_switch == True:
                        out_f.write(in_str)
                    
    return out_fpath
        
def organize_vtk(postProcPath, pattern="*.vtk"):
    print("Collecting postProcessing files into json format")

    if type(postProcPath) != pathlib.Path:
        postProcPath = pathlib.Path(postProcPath)
    resultDict = {}
    for fpath in postProcPath.rglob(pattern):
        print(fpath)
        try:
            resultDict[fpath.stem].append(dict(name=fpath.as_posix(),
                time=float(fpath.parent.name)))
        except KeyError:
            resultDict[fpath.stem] = list()
            resultDict[fpath.stem].append(dict(name=fpath.as_posix(),
                time=float(fpath.parent.name)))
    
    return resultDict

if __name__ == "__main__":

    postProcPath=pathlib.Path("./postProcessing")
    # postProcType_list = ["surfaces", "sets"]
    # out_prefix_list = ["output", "output_sets"]
    postProcType_list = ["surfaces",]
    postProcPrefix_list = ["planex0", "planey0", "planez200", "planez400",
                           "building_terrain", "two_meter_terrain",
                           "ten_meter_terrain"]
    out_prefix_list = ["output",]
    for postProcType, out_prefix in zip(postProcType_list, out_prefix_list):
        print(postProcType, out_prefix)
        postProcDest = postProcPath/postProcType
        # process each timestep into one file
        for postProcPrefix in postProcPrefix_list:
            combVtks = {}
            for t in postProcDest.iterdir():
                print(f"Combining vtk files for {t} {postProcPrefix}")
                fpaths = list(t.glob(f"{postProcPrefix}.vtp"))
                print(fpaths)
                combVtks[t.name] = combine_vtk(*fpaths,
                                               out_prefix=out_prefix+postProcPrefix)
            
            print(combVtks)
            resultDict = organize_vtk(postProcPath, pattern=postProcPrefix+".vtp")
            print(resultDict.keys())
            for key in resultDict.keys():
                print(f"Writing series file for {key}")
                with open(key + ".vtp.series", 'w') as f:
                    keyDict = {"file-series-version": "1.0", 
                                      "files" : resultDict[key]}
                    json.dump(keyDict, f)




