set -e

module load R/3.5.2 subread snakemake

read directory configfile params <<< "$@"

snakemake --snakefile Snakefile -d $directory --configfile $configfile $params
