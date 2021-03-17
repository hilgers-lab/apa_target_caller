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

* DEXseq, HTseq library (e.g. conda or pip)


## Input data

1.) PolyA database

The PolyA database needs to be a BED file, presented single nt positions of polyA sites. These polyA sites will be merged by proximity or removed if they are too close to the boundaroy of an exon segment (DEXseq jargon: exon bins). In Carrasco _et al_ 2020, we derived PolyA sites from isoSCM using RNA-seq and PARQ using 3'-seq.


* PAQR

https://github.com/zavolanlab/PAQR_KAPAC

example output files:
single_distal_sites.bed
tandem_pas_expressions.bed

* IsoSCM

https://github.com/shenkers/isoscm

example output files:
isoSCM_out.gtf




2.) Gene segmentation 

Gene segmentation happes automatically using DEXseq _dexseq_prepare_annotation.py_ which flattens the isoforms and constructs disjoin exon bins. In Carrasco _et al_ 2020, we used the Whippet segmentation, which also adds exons segment classifications, i.e. the type of splicing/APA.


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
# DEXseq path:
DEXseq_path: <DEXseq installation path>
# PolyA database path
polya_database: <polyA database>.bed
# params for breakpoint filtering
min_distance: 100
padj_cutoff: 0.05
```

Run Snakemake pipeline:

```
module load snakemake
bash run_snakemake.sh results/ config.yaml [<snakemake parameter>]
```

## Output

```
utr3_quantification/
├── Annotation
│   ├── annotation.segments.gff
│   ├── Breakpoints_pooled.merged_downstream_breakpoint.gff
│   ├── Breakpoints_pooled.merged_intervals.gff
│   ├── log
│   │   ├── Breakpoints_pooled.log
│   │   ├── exon_segmentation.log
│   │   └── Segments_split.log
│   ├── Segments_split.breakpoints_selected.gff
│   ├── Segments_split.gff
│   ├── Segments_split.nodes_selected.gff
│   └── Segments_split.saf
├── APA_targets
│   ├── <sample comparison>.APA_targets.gff
│   ├── <sample comparison>.APA_targets.locus.gff
│   ├── <sample comparison>.APA_targets.tsv
│   ├── <sample comparison>.segments_split.dexseq.gff
│   └── log
│       └── <sample comparison>.APA_targets.log
├── config.yaml
├── DEXseq
│   ├── <sample comparison>.segments_split.dexseq.tsv
│   └── log
│       └── <sample comparison>.segments_split.dexseq.log
└── featureCount
    ├── log
    │   └── utr3_quantification.segments_split.featureCounts.log
    ├── utr3_quantification.segments_split.featureCounts.tsv
    └── utr3_quantification.segments_split.featureCounts.tsv.summary


```
+ DEXSeq quanitification of each node/segment: <sample comparison>.segments_split.dexseq.tsv
+ Differential APA table: <sample comparison>.APA_targets.tsv
+ Differential APA regions: <sample comparison>.APA_targets.locus.gff
+ Segments after modification integrating PolyA database: Segments_split.gff


## Contributors

Dr. Barbara Hummel

Michael Rauer


## License
GNU GPL license (v3)
