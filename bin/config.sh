#!/bin/bash

# loaded internally
proj_file="/app/etc/proj4str.txt"
las_proj_file="/app/etc/las_proj.txt"
data_dir="/mnt/d/ubuntu_delft/workspace/chicago/DATA"
building_footprint_source="$data_dir/Buildings_20250304.csv"
las_file="$data_dir/17259075.las"
target_dir=/mnt/d/ubuntu_delft/workspace/test_result

# configuration params
lon0="-87.6372"
lat0="41.90506"
building_buffer=300
domain_buffer=400

polyprep_buffer_size="1.0"
polyprep_sim_tol="0.5"
polyprep_rm_hole="0"

# 1085000 1087499 1950000 1952499
