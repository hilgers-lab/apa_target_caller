# exaR


## Overview

[Introduction](#Introduction)

[Dependencies](#dependencies)

[Input data](#installation)

[Usage](#usage)

[Output](#output)

[Contributors](#contributors)

[License](#License)


### Introduction

exaR is a robust computational approach to quantify alternative poly(A) site usage from traditional mRNA-seq datasets.


## Dependencies

### R packages:

* readr_1.3.1
* rtracklayer_1.44.2
* GenomicFeatures
* dplyr_0.8.3
* stringr_1.4.0
* tibble_2.1.3
* ggplot2_3.2.1
* reshape2_1.4.3
* pheatmap_1.0.12
* DEXSeq_1.28.3
* janitor_2.0.1

### Other tools 

* DEXseq: 

HTseq library (e.g. conda or pip)


## Input data

1.) PAQR and/or IsoSCM

* PAQR

https://github.com/zavolanlab/PAQR_KAPAC

example output files:
single_distal_sites.bed
tandem_pas_expressions.bed

* IsoSCM

https://github.com/shenkers/isoscm

example output files:
isoSCM_out.gtf




2.) Whippet annotation

https://github.com/timbitz/Whippet.jl

example output files:
Whippet_nodes.gff


3.) sample sheet

The sample sheet is a tab-separated file with two columns, named name and condition. For comparison between two conditions, the name you assign to "condition" is not relevant, but rather the order is. The group mentioned first (in the above case "ctrl") would be used as a "control" and the group mentioned later would be used as "test".

example file:
sample_sheet.tsv


4.) gene annotation (Gencode / Ensembl / etc)

example file:
dm6_ensembl96.gtf

5.) bam files

Mapping of RNA-seq data with aligner of choice


## Usage

Include all input data paths in 'config.yaml'

```
# used as prefix
project_name: utr3_quantification
# directory with bam files
bam_dir: bam_files/
samplesheets_dir: samplesheets/
# reference annotation
annotation: dm6_ensembl96.gtf
# original whippet nodes
segments: Whippet_nodes.gff
# list of paqr/isoscm breakpoints (separated by space)
paqr:
  single_distal_sites.bed
  tandem_pas_expressions.bed
isoscm:
  isoSCM_out.gtf
# params for breakpoint filtering
min_distance: 100
isoscm_confidence: 0.7
padj_cutoff: 0.05
```

Run Snakemake pipeline:

```
module load snakemake
bash run_snakemake.sh results/ config.yaml [<snakemake parameter>]
```

## Output

```
utr3_quantification
├── Annotation
│   ├── Breakpoints_pooled.merged_downstream_breakpoint.gff
│   ├── Breakpoints_pooled.merged_intervals.gff
│   ├── Segments_split.breakpoints_selected.gff
│   ├── Segments_split.gff
│   ├── Segments_split.nodes_selected.gff
│   └── Segments_split.saf
├── APA_targets
│   ├── result.APA_targets.gff
│   ├── result.APA_targets.locus.gff
│   ├── result.APA_targets.tsv
│   ├── result.segments_split.dexseq.gff
├── config.yaml
├── DEXseq
│   ├── result.segments_split.dexseq.tsv
└── featureCount
    ├── utr3_quantification.segments_split.featureCounts.tsv
    └── utr3_quantification.segments_split.featureCounts.tsv.summary

```
+ DEXSeq quanitification of each node/segment: result.segments_split.dexseq.tsv
+ Differential APA table: result.APA_targets.tsv
+ Differential APA regions: result.APA_targets.locus.gff
+ Segments after modification using PAQR/IsoSCM: Segments_split.gff


## Contributors

Dr. Barbara Hummel

Michael Rauer


## License
GNU GPL license (v3)
