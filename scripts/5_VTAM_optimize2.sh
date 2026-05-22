#!/bin/bash

set -euo pipefail
set -x

cd $HOME/seamobb_metabarcoding_pipeline

datasets=("Arm01" "Arm02" "Arm03" "Arm04" "TEST1")

for ds in "${datasets[@]}"; do

    echo "Processing $ds"

    #################################################
    # STEP OptimizeLFNreadCountAndLFNvariant
    #################################################


    vtam optimize \
        --db db.sqlite  \
        --sortedinfo results/${ds}/2_demultiplexed/sortedinfo.tsv  \
        --sorteddir results/${ds}/2_demultiplexed \
        --known_occurrences results/${ds}/4_optimize/known_occurrences.tsv \
        --lfn_variant_replicate \
        --outdir results/${ds}/4_optimize \
        --params results/${ds}/user_input/params_optimize.yml \
        --until OptimizeLFNreadCountAndLFNvariant \
        -v 

done
