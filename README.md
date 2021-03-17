# exaR - APA target calling through exon crawling 

## Introduction

exaR is a robust computational approach to quantify alternative poly(A) site usage from traditional mRNA-seq datasets.


## Usage

Run Snakemake pipeline:

```
bash run_snakemake.sh <working directory> <config>.yaml [<snakemake parameter>]
```

Checkout the [Installation instructions](#Installation)

### Input data

The following inputs are provided by a config file

1) PolyA database: a 6-column BED file of single nucleotide coordinates defining new 3' ends of transcript isoforms, that need to be integrated into the genome annotation.
2) Genome annotation: The is the entire genome annotation (gtf) containing the full annotation, including at least gene, transcript, exon and CDS information. From this the exon segments are derived and the PolyA sites from the PolyA database are integrated.
3) Sample sheet: exaR uses DEXseq for identifying differential 3'UTR segments. Here the samplesheet needs to contains columns 'name' and 'condition', and for condition the _control is defined by the first occurrence_.
4) Alignments: These files need to be in BAM format. Can be produces using [snakePipes mRNA-seq workflow](https://snakepipes.readthedocs.io/en/latest/content/workflows/mRNA-seq.html)

### Config file:

All fields are required:

```
# used as prefix
project_name: utr3_quantification
# directory with bam files
bam_dir: bam_files/
# directory with samplessheets, tsv format and suffix
samplesheets_dir: samplesheets/
# reference annotation
annotation: dm6_ensembl96.gtf
# DEXseq path:
DEXseq_path: <DEXseq installation path>/DEXseq
# PolyA database path
polya_database: <polyA database>.bed
# params for breakpoint filtering
min_distance: 100
padj_cutoff: 0.05
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

+ DEXSeq quanitification of each node/segment: `<sample comparison>.segments_split.dexseq.tsv`
+ Differential APA table: `<sample comparison>.APA_targets.tsv`
+ Differential APA regions: `<sample comparison>.APA_targets.locus.gff`
+ Segments after modification integrating PolyA database: `Segments_split.gff`

`<sample comparison>` is the filename of the samplesheet (cropped tsv extension=


## Installation

The installation through conda can take several hours and - especially the R packages - can be installed manually as well. 

### Manual setup (Recommended)

```
conda create -n exaR  
conda activate exaR
conda install -c conda-forge r-readr r-base r-dplyr r-stringr r-tibble r-ggplot2 r-reshape2 r-pheatmap r-janitor 
conda install -c bioconda htseq snakemake subread
```

Once the conda setup is done, you can manually install the following bioconductor packages:

* [rtracklayer](https://bioconductor.org/packages/release/bioc/html/rtracklayer.html)
* [DEXseq](https://www.bioconductor.org/packages/release/bioc/html/DEXSeq.html)
* [GenomicFeatures](https://bioconductor.org/packages/release/bioc/html/GenomicFeatures.html)

Find path of _dexseq_prepare_annotation.py_

Extract libPaths from R 

```
R -e '.libPaths()'
```


and replace `<DEXseq installation path>` config entry in

```
[...]
# DEXseq path:
DEXseq_path: <DEXseq installation path>/DEXseq
[...]
```


### From conda.yaml

`conda env create -f exaR.yaml`

This conda environment config installs

* Python + HTseq
* R-base + related packages from CRAN, bioconductor *
  * Tested with R/4.0.3 
* subread for featureCounts
* snakemake 
  * Tested with python>=3.5


## Contributors

Dr. Barbara Hummel

Michael Rauer


## License

GNU GPL license (v3)
