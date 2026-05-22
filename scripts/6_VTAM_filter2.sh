#!/bin/bash

set -euo pipefail
set -x

cd $HOME/seamobb_metabarcoding_pipeline

datasets=("Arm01" "Arm02" "Arm03" "Arm04" "TEST1")

for ds in "${datasets[@]}"; do

    echo "Processing $ds"

    #################################################
    # STEP filter_default 
    #################################################
    
    mkdir -p results/${ds}/5_filter_optimized
    vtam filter \
        --db db.sqlite \
        --sortedinfo results/${ds}/2_demultiplexed/sortedinfo.tsv \
        --sorteddir results/${ds}/2_demultiplexed \
        --asvtable results/${ds}/5_filter_optimized/asvtable.tsv \
        --lfn_variant_replicate \
        --params results/${ds}/user_input/params_filter_optimized.yml \
        -v \

done
