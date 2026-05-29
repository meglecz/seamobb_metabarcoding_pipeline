library(vtamR)
library(dplyr)
library(tidyr)

setwd("~/seamobb_metabarcoding_pipeline")

###################################
# Number of sequences in input FASTQ files
rc_fastq <- count_reads_in_dir(
  dir="~/seamobb_metabarcoding_pipeline/fastq", 
  pattern="_R1_001.fastq.gz", 
  file_type="fastq"
)

sum(rc_fastq$read_count)
#[1] 51,671,544
###################################
# Number of sequences after merge, demultiplex and trim
runs <- c("Arm01", "Arm02", "Arm03", "Arm04", "TEST1" )

rc_demultiplex <- data.frame(
  filename = character(),
  read_count = numeric()
)

for(run in runs){
  
  fasta_dir <- file.path("~/seamobb_metabarcoding_pipeline/results", run, "2_demultiplexed")
  fasta_dir
  df <- count_reads_in_dir(
    dir=fasta_dir, 
    pattern="fasta.gz", 
    file_type="fasta"
  )
  
  rc_demultiplex <- rbind(rc_demultiplex, df)
}

sum(rc_demultiplex$read_count)
# [1] 37,709,707

###################################
# Number of sequences in the filtered dataset


asv_table <- read.table("results/pooled_asv_tables_ltg.tsv", sep="\t", header = TRUE, na.strings = c("NA", ""))
# delete run names before samples names
cols <- colnames(asv_table)
cols <- sub(".+\\.", "", cols)
colnames(asv_table) <- cols

rc_df <- asv_table %>%
  select(-c(marker, variant, sequence_length, read_count, sample_count, runs),-c(clusterid_7:sequence) )

sum(rc_df)
#[1] 14,419,006

###################################
# Number of sequences in the filtered dataset without external samples, and controls (negative and positive)

sample_types <- read.table("metadata/sample_types.tsv", sep="\t", header = TRUE)

asv_table_long <- asv_table %>%
  select(-c(marker, sequence_length, read_count, sample_count)) %>%
  pivot_longer(cols = -c(variant, runs, clusterid_7, clustersize, pid, ltg_taxid, ltg_name, ltg_rank, ltg_rank_index, 
                         domain, kingdom, phylum, class, order, family, genus, species, sequence ),
               names_to = "sample", 
               values_to= "read_count"
              ) %>%
  filter(read_count > 0) %>%
  left_join(sample_types, by="sample") %>%
  select(variant, sample, sample_type, mock_type, read_count, run, clusterid_7, clustersize, pid, ltg_taxid, ltg_name, 
        ltg_rank, ltg_rank_index, domain, kingdom, phylum, class, order, family, genus, species, sequence)

# delete samples from othor studies
asv_table_long <- asv_table_long %>%
  filter(!grepl(pattern="^Epi", sample ))
# Number of samples with reads
length(unique(asv_table_long$sample))
# 402

mock <- asv_table_long %>%
  filter(sample_type == "mock")
# Number of moock samples with reads
length(unique(mock$sample))
# 40


negative <- asv_table_long %>%
  filter(sample_type == "negatif")
# Number of negative controls with reads
length(unique(negative$sample))
# 5

########### dataset WO mock, negative and samples from other studies
asv_table_long <- asv_table_long %>%
  filter(sample_type != "mock" & sample_type != "negatif")
# Number of real samples with reads
length(unique(asv_table_long$sample))
# 357

########### Total number of ASV from real samples
length(unique(asv_table_long$variant))
# 6538

########### Total number of ASV from real samples assigned to at least phylum
tmp <- asv_table_long %>%
  filter(!is.na(phylum))
length(unique(tmp$variant))
# 2673

tmp <- asv_table_long %>%
  filter(!is.na(class))
length(unique(tmp$variant))
# 2306

tmp <- asv_table_long %>%
  filter(!is.na(order))
length(unique(tmp$variant))
# 1463

tmp <- asv_table_long %>%
  filter(!is.na(family))
length(unique(tmp$variant))
# 982

tmp <- asv_table_long %>%
  filter(!is.na(genus))
length(unique(tmp$variant))
# 867

tmp <- asv_table_long %>%
  filter(!is.na(species))
length(unique(tmp$variant))
# 634

################## Average number of reads per sample

tmp <- asv_table_long %>%
  group_by(sample) %>%
  summarize(read_count_sample = sum(read_count))

mean(tmp$read_count_sample)
# 30966.18
se <- sd(tmp$read_count_sample) / sqrt(length(tmp$read_count_sample))
se
# 1377.251


min(tmp$read_count_sample)
# 126
max(tmp$read_count_sample)
# 230189

# Number of FP in negative controls
FP_tneg <- nrow(negative)
FP_tneg
# 8

# Number of FP in mock

# Number of FN in mock




