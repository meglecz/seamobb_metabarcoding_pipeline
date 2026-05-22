library(dplyr)
library(yaml)

setwd("~/seamobb_metabarcoding_pipeline/")
runs = c("Arm01", "Arm02", "Arm03", "Arm04", "TEST1")

###########################################
# make directory structure

for(run in runs){
  
  dir <- file.path("results", run, "user_input")
  if(!dir.exists(dir)){
    dir.create(dir, recursive =TRUE)
  }
}

###########################################
# make fastqinfo file for each sequencing run

fastqinfo <- read.table("metadata/fastqinfo_seamobb.tsv", sep="\t", header = TRUE)
for(run_tmp in runs){
  
  outfile = file.path("results", run_tmp, "user_input", "fastqinfo.tsv")
  
  df <- fastqinfo %>%
    filter(Run == run_tmp)
  write.table(df, file = outfile, row.names = FALSE, sep="\t")
}

###########################################
# make sample_types file for each sequencing run
sample_types <- read.table("metadata/sample_types.tsv", sep="\t", header = TRUE)

for(run_tmp in runs){
  
  outfile = file.path("results", run_tmp, "user_input", "sample_types.tsv")
  
  df <- sample_types %>%
    filter(run == run_tmp)
  write.table(df, file = outfile, row.names = FALSE, sep="\t")
}

###########################################
# make mock_composition file for each sequencing run

mock_composition <- read.table("metadata/mock_composition_seamobb.tsv", sep="\t", header = TRUE)
sample_types <-sample_types %>%
  filter(run %in% runs & sample_type == "mock") %>%
  select(run, sample, mock_type)

for(run_tmp in runs){

  outfile = file.path("results", run_tmp, "user_input", "mock_composition.tsv")
  
  mock_tmp <- sample_types %>%
    filter(run == run_tmp) %>%
    left_join(mock_composition, by="mock_type", relationship="many-to-many") %>%
    mutate( mock=1, variant=NA, tax_name=taxon) %>%
    select(marker, run, sample, mock, variant, action, sequence, tax_name)

  write.table(mock_tmp, file = outfile, row.names = FALSE, sep="\t")
}

###########################################
# yml parameter files for optimize2 for each sequencing run

options(scipen = 999)
params <- read.table("metadata/parameters.tsv", sep="\t", header = TRUE)
for(run_tmp in runs){
  
  outfile = file.path("results", run_tmp, "user_input", "params_optimize.yml")
  
  tmp <- params %>%
    filter(run == run_tmp) %>%
    select(global_read_count_cutoff, pcr_error_var_prop, lfn_sample_replicate_cutoff)
  
  tmp <- t(as.matrix(tmp))
  

  params_list <- as.list(tmp[,1])
  names(params_list) <- rownames(tmp)
  
  # convert selected parameter to integer
  params_list$global_read_count_cutoff <- as.integer(
    params_list$global_read_count_cutoff
  )

  write_yaml(params_list, file = outfile)
}

###########################################
# yml parameter files for filter2 for each sequencing run

for(run_tmp in runs){
  
  outfile = file.path("results", run_tmp, "user_input", "params_filter_optimized.yml")
  
  tmp <- params %>%
    filter(run == run_tmp) %>%
    select(-run)
  
  tmp <- t(as.matrix(tmp))
  
  params_list <- as.list(tmp[,1])
  names(params_list) <- rownames(tmp)
  
  # convert selected parameter to integer
  params_list$global_read_count_cutoff <- as.integer(
    params_list$global_read_count_cutoff
  )
  params_list$lfn_read_count_cutoff <- as.integer(
    params_list$lfn_read_count_cutoff
  )

  write_yaml(params_list, file = outfile)
}


