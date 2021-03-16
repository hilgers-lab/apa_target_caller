set -e

module load R/4.0.3 subread snakemake

read directory configfile params <<< "$@"

snakemake --snakefile Snakefile -d $directory --configfile $configfile $params
