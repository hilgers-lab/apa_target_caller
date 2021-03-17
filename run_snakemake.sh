set -e

module load subread R/4.0.3

read directory configfile params <<< "$@"

snakemake --snakefile Snakefile -d $directory --configfile $configfile $params
