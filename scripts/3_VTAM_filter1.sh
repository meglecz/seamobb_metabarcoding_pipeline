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
    
    mkdir -p results/${ds}/3_filter_default
    vtam filter \
        --db results/db.sqlite \
        --sortedinfo results/${ds}/2_demultiplexed/sortedinfo.tsv \
        --sorteddir results/${ds}/2_demultiplexed \
        --asvtable results/${ds}/3_filter_default/asvtable_default.tsv \
        --lfn_variant_replicate \
        --params metadata/params_filter_default.yml \
        -v \

done
