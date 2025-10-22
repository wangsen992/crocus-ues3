#!/bin/bash

# # loaded internally
# proj_file="/app/etc/proj4str.txt"
# data_dir="/mnt/d/ubunt_delft/workspace/chicago/DATA"
# building_footprint_source="$data_dir/Buildings_20250304.csv"
# las_file="$data_dir/17259075.las"
# target_dir=/mnt/d/ubuntu_delft/workspace/test_result
# 
# # configuration params
# lon0="-87.6372"
# lat0="41.90506"
# building_buffer=100
# domain_buffer=300

source $1

echo $NSLOTS

the_source=$(readlink -f -- "${BASH_SOURCE[0]}")
the_dirname=$(dirname "${the_source}")

# create proj file with the pre-set lat lon
echo "+proj=lcc +lat_0=$lat0 +lon_0=$lon0 +lat_1=33 +lat_2=45 +ellps=GRS80" >> $target_dir/proj4str.txt
proj_file=$target_dir/proj4str.txt

# run create_bnbs
echo "Creating building and domain boundary geojson files..."
python3 $the_dirname/create_bnds.py \
        --proj_fname $proj_file \
        --lon0 $lon0 \
        --lat0 $lat0 \
        --building_buffer $building_buffer \
        --domain_buffer $domain_buffer \
        --target_dir $target_dir

# Run directly after create_bnbs to use the influenceRegion.json
echo "Subsetting building footprint files from source file..."
python3 $the_dirname/subset_bldfprt.py \
        --proj_fname $proj_file \
        --source_bldfprt $building_footprint_source \
        --target_dir $target_dir

# Run polyprep to simplify 
echo "Simplifying buildings geojson..."
python3 /app/City4CFD/tools/polyprep/polyprep.py $target_dir/buildings.geojson $target_dir/buildings_simplified.geojson $polyprep_buffer_size --simplification_tol $polyprep_sim_tol --remove_holes $polyprep_rm_hole

# subset las files to domain region
echo "Subsetting big las file to smaller size"
xy_region=$(python3 $the_dirname/get_las_bnbs.py \
										--subset_geojson $target_dir/domainBnd.geojson \
										--proj_fname $proj_file \
										--las_proj_fname $las_proj_file )
subset_las_file=$target_dir/subset_lidar.las
echo $xy_region
mkdir $target_dir/subset_las
names=$(ls $data_dir/lidar_cook/cook-las2/*.las | xargs -I {} basename {})
echo $names | \
	xargs -d " " -P 4 -I {} /app/LAStools/bin/las2las64 \
		-i $data_dir/lidar_cook/cook-las2/{} -o $target_dir/subset_las/subset_{} \
		-inside $xy_region
# /app/LAStools/bin/lasmerge64 -i $target_dir/subset_las/*.las -o $target_dir/subset_lidar.las
 
# /app/LAStools/bin/las2las64 -i $las_file -o $subset_las_file -inside $xy_region
    
# create las files
echo "Preprocessing source las files to target crs and separate into building and ground..." 
python3 $the_dirname/prep_las_par.py \
    --proj_fname $proj_file \
    --source_las $target_dir/subset_las \
    --subset_geojson $target_dir/influenceRegion.geojson \
    --target_dir $target_dir \
		--num_workers $NSLOTS

# Now merge all las files accordingly
echo "merge subset las files in each category"
/app/LAStools/bin/lasmerge64 -i $target_dir/ground_las/*.las -o $target_dir/ground_raw.las &
/app/LAStools/bin/lasmerge64 -i $target_dir/building_las/*.las -o $target_dir/building_raw.las &
/app/LAStools/bin/lasmerge64 -i $target_dir/vegetation_las/*.las -o $target_dir/vegetation_raw.las &
/app/LAStools/bin/lasmerge64 -i $target_dir/water_las/*.las -o $target_dir/water.las &
wait

echo "Subsample las files in each category"
/app/LAStools/bin/las2las64 -i $target_dir/ground_raw.las -o $target_dir/ground.las -keep_random_fraction 0.05 &
/app/LAStools/bin/las2las64 -i $target_dir/building_raw.las -o $target_dir/building.las -keep_random_fraction 0.05 &
/app/LAStools/bin/las2las64 -i $target_dir/vegetation_raw.las -o $target_dir/vegetation.las -keep_random_fraction 0.05 &
wait

echo "Complete run.sh file" 
echo
echo
