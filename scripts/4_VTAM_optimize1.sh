#!/bin/bash

set -euo pipefail
set -x

cd $HOME/seamobb_metabarcoding_pipeline
#conda activate vtam_pipeline

datasets=("Arm01" "Arm02" "Arm03" "Arm04" "TEST1")

for ds in "${datasets[@]}"; do

    echo "Processing $ds"

    #################################################
    # STEP make_known_occurrences
    #################################################
    
    mkdir -p results/${ds}/4_optimize
    vtam make_known_occurrences \
        --asvtable  results/${ds}/3_filter_default/asvtable_default.tsv \
        --sample_types results/${ds}/user_input/sample_types.tsv \
        --mock_composition results/${ds}/user_input/mock_composition.tsv \
        --known_occurrences results/${ds}/4_optimize/known_occurrences.tsv \
        --missing_occurrences results/${ds}/4_optimize/missing_occurrences.tsv \
        -v

    #################################################
    # STEP OptimizePCRerror 
    #################################################

    vtam optimize \
        --db db.sqlite  \
        --sortedinfo results/${ds}/2_demultiplexed/sortedinfo.tsv  \
        --sorteddir results/${ds}/2_demultiplexed \
        --known_occurrences results/${ds}/4_optimize/known_occurrences.tsv \
        --lfn_variant_replicate \
        --outdir results/${ds}/4_optimize \
        --until OptimizePCRerror \
        -v 
    
    #################################################
    # STEP OptimizeLFNsampleReplicate
    #################################################

    vtam optimize \
        --db db.sqlite  \
        --sortedinfo results/${ds}/2_demultiplexed/sortedinfo.tsv  \
        --sorteddir results/${ds}/2_demultiplexed \
        --known_occurrences results/${ds}/4_optimize/known_occurrences.tsv \
        --lfn_variant_replicate \
        --outdir results/${ds}/4_optimize \
        --until OptimizeLFNsampleReplicate \
        -v 
done
