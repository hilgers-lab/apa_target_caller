suppressPackageStartupMessages(library(optparse))

# specify our desired options in a list
# by default OptionParser will add an help option equivalent to 
# make_option(c("-h", "--help"), action="store_true", default=FALSE, 
#               help="Show this help message and exit")
option_list <- list( 
  # make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
  #             help="Print extra output [default]"),
  # make_option(c("-q", "--quietly"), action="store_false", 
  #             dest="verbose", help="Print little output"),
  make_option(c("-g", "--gtf"), type="character",
              help="Reference annotation (ensembl gtf)"),
  make_option(c("-s", "--segments"), type="character", 
              help = "Segments, i.e. Whippet nodes (gtf/gff format)"),
  make_option(c("-b", "--breakpoints"), type="character", 
              help = "Breakpoint predictions (gtf format)"),
  
  make_option(c("-m", "--method"), type="character", default = 'DEXseq',
              help = "Method used to create segments. \'DEXseq\' or \'Whippet\' (default: %default)"),
  
  make_option("--min.distance", type="integer",  default = 100, 
              help = "Minimum distance of breakpoint to the node's edges (default: %default)"),
  
  make_option(c("-d", "--debug"), type="logical", action = 'strore_true', default = FALSE, 
              help = "Output intermediate files (default: %default)"),
  
  make_option(c("-o", "--outfileName"), type="character", 
              help = "Name of file output file (gff). Accompanied by an saf file (for featureCounts)")
)

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
opt <- parse_args(OptionParser(option_list=option_list))
opt$prefix = gsub('.gff$', '',opt$outfileName)

assertthat::assert_that(opt$method %in% c('DEXseq','Whippet'))

suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))


addDistanceInNode <- function(breakpoints, nodes) {
  require(GenomicFeatures)
  d.end <- distanceToNearest(breakpoints, resize(nodes, fix = 'end', width = 1))
  d.start <- distanceToNearest(breakpoints, resize(nodes, fix = 'start', width = 1))
  
  d.min <- apply(cbind(mcols(d.end)$distance, mcols(d.start)$distance), 1, min)
  
  mcols(breakpoints)$distance = d.min
  return(breakpoints)
}


params = list(min.dist = opt$min.distance)

# utr3 as reference
txdb0 <- makeTxDbFromGFF(opt$gtf)
dm6.utr3 <- threeUTRsByTranscript(txdb0, use.names = TRUE)
dm6.tes <-  resize(transcripts(txdb0, use.names = TRUE), fix = 'end', width = 1)

dm6.ens96 <- import.gff(opt$gtf, feature.type = 'gene')
biotype.map <- deframe(mcols(dm6.ens96)[c('gene_id','gene_biotype')])

# Import and filter whippet nodes
breakpoints_all <- import.gff(opt$breakpoints)
segments <- import.gff(opt$segments)

if(opt$method == 'DEXseq'){
  mcols(segments)$Gene = mcols(segments)$gene_id
  mcols(segments)$Node = mcols(segments)$exonic_part_number
}
# if(opt$method == 'Whippet'){
#   # segments
# }



params$chrom.valid = as.character(unique(intersect(seqnames(segments), seqnames(breakpoints_all))))

seqlevels(breakpoints_all, pruning.mode = 'coarse') <- params$chrom.valid
seqlevels(segments, pruning.mode = 'coarse') <- params$chrom.valid


mcols(segments)$is.utr3 <- overlapsAny(segments, dm6.utr3)
segments_nodes <- subset(segments, 
                         biotype.map[Gene] %in% c('protein_coding','pseudogene') & 
                           seqnames %in% params$chrom.valid)


# Import and filter breakpoints
## Filter A) by: 
# - min distance
# - utr3 overlapping 
# - far enough from nearest TES
dd.tes <- distanceToNearest(breakpoints_all, dm6.tes)
mcols(breakpoints_all)$distance2tes <- Inf
mcols(breakpoints_all)$distance2tes[queryHits(dd.tes)] <- mcols(dd.tes)$distance

breakpoints <- subset(breakpoints_all, 
                      overlapsAny(breakpoints_all, dm6.utr3) & 
                        distance2tes > params$min.dist)

if(length(breakpoints) == 0)
  stop('No breakpoints left after filtering.\nExit...')

# Filter B)
# require min distance from node start/end
breakpoints <- addDistanceInNode(breakpoints, segments_nodes)

breakpoints.filtered <- subset(breakpoints, distance > params$min.dist)

ii <- findOverlaps(breakpoints.filtered, segments_nodes)
# Note: missed breakpoints are from ambiguous whippet nodes

# for each node, cut to pieces by breakpoints
cutNode <- function(start, end, node.template){
  node.template$width = NULL
  node.template$start = start
  node.template$end = end
  return(node.template)
}

cat('Breaking nodes\n')
segment.mat <- as.data.frame(segments_nodes)
segment_split <- data.frame()
for(j in sort(unique(subjectHits(ii)))){
  # print(j)
  node.x <- segment.mat[j,]
  i <- queryHits(subset(ii, subjectHits == j))
  breaks.start <- c(node.x$start, start(breakpoints.filtered[i]))
  breaks.end <- c(breaks.start[-1], node.x$end)
  
  node.y <- node.x
  nodes_new <- do.call(rbind, mapply(cutNode, breaks.start, breaks.end, MoreArgs = list(node.template=node.x), SIMPLIFY = FALSE))
  segment_split <- rbind(segment_split, nodes_new)
}

cat('Negative width nodes:\n')
segment_split[which(segment_split$end - segment_split$start < 0),]


segments_tmp <- makeGRangesFromDataFrame(segment_split, keep.extra.columns = TRUE)
tmp = as.data.frame(segments[!segments %over% segments_tmp])
tmp$width = NULL
segments_novel <- sort(makeGRangesFromDataFrame(rbind(tmp,
                                                      segment_split), keep.extra.columns = TRUE))

tab0 <- as.data.frame(segments_novel)
tab1 <- tab0 %>% 
  group_by(Gene) %>% 
  add_tally() %>%
  mutate(Node = ifelse(strand == '+', order(start), n - order(start) + 1))

tab.out <- data.frame(name = paste(tab1$Gene, tab1$Node, tab1$Type, sep = '_'), tab1)
assertthat::assert_that(all(table(tab.out$name) == 1 ), msg = 'Something went wrong. There are nodes in a gene with the same index. \nExit...')

# nodes to breaks
if(opt$debug){
  segments2break <- segments_nodes[unique(subjectHits(ii))]
  export.gff(segments2break,paste0(opt$prefix,'.nodes_selected.gff'))
  # selected breakpoints
  breakpoints.selected <- breakpoints.filtered[queryHits(ii)]
  export.gff(breakpoints.selected, paste0(opt$prefix,'.breakpoints_selected.gff'))
}

write_tsv(tab.out, paste0(opt$prefix,'.saf'), col_names = FALSE)

segments_novel <- makeGRangesFromDataFrame(tab.out, keep.extra.columns = TRUE)
export.gff3(segments_novel, paste0(opt$prefix,'.gff'))

