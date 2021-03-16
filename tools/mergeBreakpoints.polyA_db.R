library(optparse)

option_list <- list( 
  make_option(c("-g", "--gtf"), type="character", 
              help="Reference annotation e.g ensembl annotation (gtf)"),
  make_option(c("-s", "--segments"), type="character", 
              help="Exon segments, e.g. Whippet node output (gtf)"),
  make_option(c("-p", "--polya_db"), type="character",
              help = "polyA database in BED format"),
  
  make_option("--minimum.distance", type="integer", default = 100,
              help = "Minimum distance between breakpoints and breakpoint to node edge (default: %default)"),
  make_option("--isoscm.confidence", type="double", default = 0.7,
              help = "Confidence cutoff for isoscm nodes (default >= %default)"),
  
  make_option(c("-o", "--outfileNamePrefix"), type="character",
              help = "Prefix for output gtf. Files are '<prefix>.merged_intervals.gff' and '<prefix>.merged_downstream_breakpoint.gff'")
)

opt <- parse_args(OptionParser(option_list=option_list))
print(opt)

suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(GenomicFeatures))

chrom.filter <- function(gr, chroms.valid){
  require(GenomicFeatures)
  seqlevels(gr, pruning.mode = 'coarse') = chroms.valid
  return(gr)
}

addDistanceInNode <- function(breakpoints, nodes) {
  require(GenomicFeatures)
  d.end <- distanceToNearest(breakpoints, resize(nodes, fix = 'end', width = 1))
  d.start <- distanceToNearest(breakpoints, resize(nodes, fix = 'start', width = 1))
  
  d.min <- apply(cbind(mcols(d.end)$distance, mcols(d.start)$distance), 1, min)
  
  mcols(breakpoints)$distance = d.min
  return(breakpoints)
}


import.polya_db <- function(paqr_file){
  return(import.bed(paqr_file))
}


filter.chromosomes <- function(gr, chrom.valid){
  return(chrom.filter(gr, chrom.valid))
}

filter.breakPoints <- function(breakpoints, segments, min.distance = 100){
  breakpoints = addDistanceInNode(breakpoints, nodes = segments)
  return(subset(breakpoints, distance >= min.distance))
}


params = list(min.distance = opt$minimum.distance,
              max.gap = opt$minimum.distance,
              min.confidence = opt$isoscm.confidence)

cat("Importing data:\n")

txdb0 = makeTxDbFromGFF(opt$gtf)
annotation.utr3 = threeUTRsByTranscript(txdb0, use.names = TRUE)

segments = import.gff(opt$segments)

cat(">> Import PolyA Database:\n")
breakpoints <- import.polya_db(opt$polya_db)


# valid chromosomes: common to all external sets - and chromsomes that contain genes with UTR3
chroms.utr3 <- as.character(subset(as.data.frame(table(seqnames(unlist(annotation.utr3)))), Freq > 0)$Var)
chrom.intersect <- intersect(intersect(unique(unlist(seqlevels(breakpoints))), seqlevels(segments)), chroms.utr3)

# filter breakpoints
breakpoints = filter.chromosomes(breakpoints, chrom.intersect)
segments = filter.chromosomes(segments, chrom.intersect)

breakpoints = filter.breakPoints(breakpoints , segments, params$min.distance)


# merge breakpoints
breakpoints_pooled = granges(breakpoints)
breakpoints_merged = reduce(breakpoints_pooled, min.gapwidth = params$max.gap)

mcols(breakpoints_merged)$is.clustered = width(breakpoints_merged) > 1
cat('---\n')
cat('- Number of total breakpoints: ', length(breakpoints_pooled),'\n')
cat('- Number of merge breakpoints: ', length(breakpoints_merged),'\n')
cat('- Number of clustered breakpoints: ', sum(mcols(breakpoints_merged)$is.clustered),'\n')
cat('---\n')

breakpoints.final = resize(breakpoints_merged, fix = 'end', width = 1)


# export breakpoints.final
cat("Writing files to\n")
cat(opt$outfileNamePrefix, '\n')
cat(' -- ', paste0(opt$outfileNamePrefix,'.merged_intervals.gff'),'\n')
export.gff(breakpoints_merged, paste0(opt$outfileNamePrefix,'.merged_intervals.gff'))
cat(' -- ',paste0(opt$outfileNamePrefix,'.merged_downstream_breakpoint.gff'), '\n')
export.gff(breakpoints.final, paste0(opt$outfileNamePrefix,'.merged_downstream_breakpoint.gff'))
