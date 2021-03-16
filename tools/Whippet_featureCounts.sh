module load subread

read outfile annotation bams <<< "$@"

nproc=16
# bams="$(ls bams/*.bam | paste -sd ' ')"
# annotation="Whippet_nodes.saf"

# outfile="Whippet_nodes.featureCounts.tsv"

params="-t sequence_feature -g feature_id -f -s 2 --primary -p -B -T $nproc"
featureCounts $params -a $annotation -F SAF -o $outfile $bams
