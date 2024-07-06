#!/bin/bash
iwd=$(pwd) # iwd is the initial working directory
raw_data_dir = $(pwd)/raw_data
sdf_dir="$(pwd)/mols"
fp_dir="$(pwd)/fps"
mkdir -p "$sdf_dir" # will not create if already exists
mkdir -p "$fp_dir" 
script_dir="$(dirname "$0")/0_input_preparation"

# for each type of moldata the executions will be in parrallel, where the fps will be calculated after successful sdf generation
for script in "$script_dir"/[!get]*.py; do 
    script_name=$(basename -s .py "$script")
    sdf_file="$sdf_dir/$(basename -s .py "$script").sdf"
    cd "$sdf_dir" && python "$script" "$raw_data_dir" >> sdf_prep.log && cd "$fp_dir" && python "get_fps.py" "$sdf_file" >> fp_prep.log & 
done

wait
cd "$iwd"