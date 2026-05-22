# Bioinformatics Protocol

Protocol associated with the paper:

**Separating faces in ARMS metabarcoding improves marine biodiversity monitoring: a comparison across protocols, experimental designs, and photographic surveys**
Chenuil *et al*., 2026 (DOI: xxxxxxxxxxxx)

---


# Requirements

## Software

The pipeline relies on the following software and tools:

* **R (≥ 4.0)** with packages  `yaml` and `dplyr`
  Used for metadata preparation and parameter file generation.

* **VTAM (0.2.0)**
  Metabarcoding analysis pipeline used for merging, demultiplexing, filtering, and ASV inference.
  Documentation: [https://vtam.readthedocs.io/en/latest/index.html](https://vtam.readthedocs.io/en/latest/index.html) 

* **mkCOInr (≥ v0.1)**
  Toolkit used to build and customize the COInr reference database.
  Documentation: [https://mkcoinr.readthedocs.io/en/latest/](https://mkcoinr.readthedocs.io/en/latest/)

---

## Hardware requirements (recommended)

* ≥ 16 GB RAM (32 GB recommended for full runs)
* Multi-core CPU (≥ 4 cores recommended)
* ≥ 50 GB free disk space (depending on raw data size)

---

# Data

Each Illumina MiSeq sequencing run contains 96 samples, including environmental samples, negative controls, and mock communities. Each sample is amplified by PCR in triplicate. One replicate set of the 96 samples forms a sequencing library identified by MIDs. Therefore, a sequencing run produces six FASTQ files: forward and reverse reads for each of the three PCR replicates.

Each pair of FASTQ files must be demultiplexed according to the sequencing tags added to the amplicons during PCR, and the tags and primers must then be trimmed.

VTAM first merges paired-end FASTQ sequences and converts them into FASTA files, which are subsequently demultiplexed and trimmed of tags and primers.

This procedure is not directly compatible with SRA archiving, since SRA requires paired FASTQ files that are already demultiplexed and trimmed. To enable submission to SRA, demultiplexing and trimming were performed before submission; however, this preprocessing step is independent from the VTAM analyses described below.

---

## FASTQ files

The pipeline below describes the analyses performed starting from non-demultiplexed FASTQ files. It is therefore not directly applicable to the demultiplexed SRA files.

The raw non-demultiplexed FASTQ files used as input for the pipeline can be downloaded from Zenodo: [10.5281/zenodo.20344491](10.5281/zenodo.20344491)


---

## Metadata

### `fastqinfo_seamobb.tsv`

Contains information on tag combinations used for sample demultiplexing.

Columns:

* `TagFwd` — forward tag sequence
* `PrimerFwd` — forward primer sequence
* `TagRev` — reverse tag sequence
* `PrimerRev` — reverse primer sequence
* `Marker` — sequencing marker corresponding to the amplified COI region
* `Sample` — sample name
* `Replicate` — PCR replicate ID
* `Run` — sequencing run
* `FastqFwd` — forward FASTQ file name
* `FastqRev` — reverse FASTQ file name

---

### `mock_composition_seamobb.tsv`

Contains the expected ASVs for the species composing each mock community.

Columns:

* `mock_type` — name of the mock community (multiple samples may share the same mock type)
* `taxon` — taxon name
* `marker` — marker name
* `action` — either:

  * `keep`: sequences expected and retained in the dataset
  * `tolerate`: sequences that may occur accidentally at low abundance and not necessarily amplify in every mock sample
* `sequence` — expected ASV sequence

---

### `sample_types.tsv`

List of all samples.

Columns:

* `run` — sequencing run
* `sample` — sample name
* `sample_type` — `real`, `mock`, or `negative`
* `habitat` — habitat type used to identify potential false positives (e.g. ALI mock taxa are terrestrial organisms and are therefore not expected in marine environmental samples)
* `mock_type` — name of the mock community
* `marker` — sequencing marker corresponding to the amplified COI region

---

### `parameters.tsv`

Contains parameter values used in the YAML configuration files for the `optimize2` and `filter2` VTAM steps.

These values are provided to facilitate reproducibility of the pipeline.

* `pcr_error_var_prop` and `lfn_sample_replicate_cutoff` were determined after the `optimize1` step.
* `lfn_read_count_cutoff` and `lfn_variant_replicate_cutoff` were determined after the `optimize2` step.

---

## Shared parameter files

The following parameter files are identical for all sequencing runs:

* `params_filter_default.yml`
* `params_merge.yml`
* `params_sortreads.yml`

---

## Additional files

### `taxon_list_insecta.tsv`

List of taxa removed from the COInr database during database customization.

### `private_sequences.tsv`

COI reference sequences generated in the laboratory and unpublished as of 2022.

### `pool_runs.tsv`

TSV file listing the run-marker combinations to pool after filtering.

---

# Reference Database

A custom reference database for taxonomic assignment was generated from the COInr database [version 2022-05-06](https://zenodo.org/records/6555985) using [mkCOInr](https://zenodo.org/records/6566165) v0.1 ([Meglécz, 2023](https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.13756)).

COInr contains COI sequences from both BOLD and GenBank and is dereplicated within taxa.

Insect sequences were removed because they represent a major proportion of COInr but are not relevant for marine biodiversity studies. Additional unpublished reference sequences generated in the laboratory were also added.

The custom database can be downloaded from Zenodo:

[10.5281/zenodo.20344491](10.5281/zenodo.20344491)

The commands below describe the database construction process.

---

## Download COInr_2022_05_06

Download `COInr_2022_05_06.tar.gz` from:

[https://zenodo.org/records/6555985](https://zenodo.org/records/6555985)

Then extract the archive:

```bash
cd ~/seamobb_metabarcoding_pipeline
tar -zxvf COInr_2022_05_06.tar.gz
```

---

## Remove insect sequences

```bash
perl ~/mkCOInr/scripts/select_taxa.pl \
    -taxon_list metadata/taxon_list_insecta.tsv \
    -tsv COInr_2022_05_06/COInr.tsv \
    -taxonomy COInr_2022_05_06/taxonomy.tsv \
    -outdir COInr_2022_05_06_seamobb/1_WO_insecta \
    -out COInr_WO_insecta.tsv \
    -negative_list 1
```

---

## Add private and mock sequences

### Format custom sequences

```bash
perl ~/mkCOInr/scripts/format_custom.pl \
    -custom metadata/private_sequences.tsv \
    -taxonomy COInr_2022_05_06/taxonomy.tsv \
    -outdir COInr_2022_05_06_seamobb/2_custom/1_format
```

Check and manually curate the `custom_lineages.tsv` file:

* If the `homonymy` column contains `1`, remove lines corresponding to incorrect taxa.
* Verify that the lineages suggested by mkCOInr are biologically coherent.
* Save the curated file as `custom_lineages_verified.tsv`.

---

### Add taxids to custom sequences

```bash
perl ~/mkCOInr/scripts/add_taxids.pl \
    -lineages COInr_2022_05_06_seamobb/2_custom/1_format/custom_lineages_verified.tsv \
    -sequences COInr_2022_05_06_seamobb/2_custom/1_format/custom_sequences.tsv \
    -taxonomy COInr_2022_05_06/taxonomy.tsv \
    -outdir COInr_2022_05_06_seamobb/2_custom/2_add_taxids
```

---

### Move updated taxonomy

```bash
mv COInr_2022_05_06_seamobb/2_custom/2_add_taxids/taxonomy_updated.tsv \
   COInr_2022_05_06_seamobb/2_custom/taxonomy_COInr_seamobb_2022.tsv
```

---

### Dereplicate custom sequences

```bash
perl ~/mkCOInr/scripts/dereplicate.pl \
    -tsv COInr_2022_05_06_seamobb/2_custom/2_add_taxids/sequences_with_taxIDs.tsv \
    -outdir COInr_2022_05_06_seamobb/2_custom/3_dereplicate \
    -out custom_dereplicated_sequences.tsv
```

---

### Pool and dereplicate databases

```bash
perl ~/mkCOInr/scripts/pool_and_dereplicate.pl \
    -tsv1 COInr_2022_05_06_seamobb/1_WO_insecta/COInr_WO_insecta.tsv \
    -tsv2 COInr_2022_05_06_seamobb/2_custom/3_dereplicate/custom_dereplicated_sequences.tsv \
    -outdir COInr_2022_05_06_seamobb/2_custom/ \
    -out COInr_seamobb_2022.tsv
```

---

## Format database for VTAM

```bash
perl ~/mkCOInr/scripts/format_db.pl \
    -tsv COInr_2022_05_06_seamobb/2_custom/COInr_seamobb_2022.tsv \
    -taxonomy COInr_2022_05_06_seamobb/2_custom/taxonomy_COInr_seamobb_2022.tsv \
    -outfmt vtam \
    -outdir COInr_2022_05_06_seamobb/COInr_seamobb_2022 \
    -out COInr_seamobb_2022
```

---

## Clean up

```bash
mv COInr_2022_05_06_seamobb/COInr_seamobb_2022 COInr_seamobb_2022

rm -rI COInr_2022_05_06_seamobb
rm -rI COInr_2022_05_06
```

---

# VTAM Pipeline

Sequencing runs are analysed independently and pooled after filtering.

VTAM is run in a conda environnement:

~~~
conda activate vtam
~~~



---

## Prepare metadata files for each run

The script `1_make_individual_metafiles.R` generates individual:

* `fastqinfo`
* `mock_composition`
* `sample_types`
* `params_optimize.yml`
* `params_filter_optimized.yml`

files for each sequencing run.

```bash
cd ~/seamobb_metabarcoding_pipeline

Rscript scripts/1_make_individual_metafiles.R
```

---

## Merge and demultiplex

For each sequencing run:

1. Paired-end FASTQ reads are merged and converted to FASTA format using the VTAM `merge` command.
2. Merged FASTA files are demultiplexed and trimmed using the VTAM `sortreads` command.

These steps can be run for all sequencing runs using:

```bash
bash scripts/2_VTAM_merge_sortreads.sh
```

---

## Prefiltering

Initial VTAM filtering is performed using default parameters.

```bash
bash scripts/3_VTAM_filter1.sh
```

---

## Parameter optimization

The VTAM `make_known_occurrences` command identifies:

* true positives
* false positives
* false negatives

based on mock community composition negative controls and the occurrence of ASVs in incompatible habitats.

These occurrences are used to optimize the following parameters:

* `lfn_sample_replicate_cutoff`
* `pcr_error_var_prop`
* `lfn_variant_replicate_cutoff`
* `lfn_read_count_cutoff`

---

### Optimize PCR error and sample replicate cutoffs

The values of `lfn_sample_replicate_cutoff` and `pcr_error_var_prop` are first optimized for each sequencing run.

```bash
bash scripts/4_VTAM_optimize1.sh
```

The outputs should be examined manually and the parameter values validated by the user before generating the corresponding `params_optimize.yml` files.

For reproducibility and automation, these files are automatically generated by the `1_make_individual_metafiles.R` script. However, for real-world analyses, manual inspection and adjustment are recommended.

---

### Optimize low-frequency noise parameters

The optimized values above are then used to determine the best combination of:

* `lfn_variant_replicate_cutoff`
* `lfn_read_count_cutoff`

```bash
bash scripts/5_VTAM_optimize2.sh
```

The outputs should again be inspected manually before generating the final `params_filter_optimized.yml` files.

As above, these files are automatically generated for reproducibility purposes but would normally require manual curation in standard analyses.

---

## Filtering with optimized parameters

Filtering is rerun using the optimized parameter values.

```bash
bash scripts/6_VTAM_filter2.sh
```

---

## Pool sequencing runs and assign taxonomy

The filtered results from the five sequencing runs are pooled into a final ASV table, and ASVs are taxonomically assigned.

```bash
bash scripts/7_VTAM_pool_taxassign.sh
```



