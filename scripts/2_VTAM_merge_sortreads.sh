#!/bin/bash

set -euo pipefail
set -x

cd $HOME/seamobb_metabarcoding_pipeline
datasets=("Arm01" "Arm02" "Arm03" "Arm04" "TEST1")

for ds in "${datasets[@]}"; do

    echo "Processing $ds"

    #################################################
    # STEP merge
    #################################################
    
    mkdir -p results/${ds}/1_merge
    vtam merge \
        --fastqinfo results/${ds}/user_input/fastqinfo.tsv \
        --fastqdir fastq \
        --fastainfo results/${ds}/1_merge/fastainfo.tsv \
        --fastadir results/${ds}/1_merge \
        --params metadata/params_merge.yml \
        -v 


    #################################################
    # STEP sortreads 
    #################################################

    mkdir -p results/${ds}/2_demultiplexed
    vtam sortreads \
        --fastainfo results/${ds}/1_merge/fastainfo.tsv \
        --fastadir results/${ds}/1_merge \
        --sorteddir results/${ds}/2_demultiplexed \
        --params metadata/params_sortreads.yml \
        -v 

done
