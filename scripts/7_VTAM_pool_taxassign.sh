#!/bin/bash

set -euo pipefail
set -x

cd $HOME/seamobb_metabarcoding_pipeline

db_dir="COInr_seamobb_2022"
db_name="COInr_seamobb_2022"
taxonomy="COInr_seamobb_2022/COInr_seamobb_2022_taxonomy.tsv"


#################################################
# STEP Pool the filtered results of all runs
#################################################

vtam pool \
    --db db.sqlite \
    --runmarker metadata/pool_runs.tsv \
    --asvtable results/pooled_asvtable.tsv \
    -v

#################################################
# STEP taxassign
#################################################

vtam taxassign \
    --db db.sqlite \
    --asvtable results/pooled_asvtable.tsv \
    --output results/pooled_asvtable_taxa.tsv \
    --taxonomy $taxonomy \
    --blastdbdir $db_dir \
    --blastdbname $db_name \
    -v 



