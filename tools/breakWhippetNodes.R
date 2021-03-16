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
  
  make_option("--isoscm.confidence", type="double", default = 0.7,
              help = "Minimum for isoSCM's confidence (default: %default)"),
  make_option("--min.distance", type="integer",  default = 100, 
              help = "Minimum distance of breakpoint to the node's edges (default: %default)"),
  
  make_option(c("-d", "--debug"), type="logical", action = 'strore_true', default = FALSE, 
              help = "Output intermediate files (default: %default)"),
  
  make_option(c("-o", "--outfileName"), type="character", 
              help = "Name of file output file (gff). Accompanied by an saf file (for featureCounts)")
)

# nodes to breaks
# export.gff(whippet2break,paste0(opt$prefix,'.nodes_selected.gff'))
# # selected breakpoints
# export.gff(breakpoints.selected,paste0(opt$prefix,'.breakpoints_selected.gff'))
# write_tsv(tab.out, paste0(opt$prefix,'.saf'), col_names = FALSE)
# export.gff3(whippet_novel, paste0(opt$prefix,'.gff'))
# 


# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
opt <- parse_args(OptionParser(option_list=option_list))
opt$prefix = gsub('.gff$', '',opt$outfileName)

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


### Test data
# opt$gtf="data/dm6_ensembl96.gtf"
# opt$segments="data/Whippet_nodes.gff"
# opt$breakpoints="data/compare_elav_pooled_all_wt_pooled_all.gtf"
# opt$isoscm.confidence=0.7
# opt$min.distance=300

# opt$prefix="data/dm6_ensembl96.novel_breakspoints"


# opt$gtf = '/data/hilgers/group/Paper_ELAV-dependentTargets/Annotations/dm6_ensembl96.gtf'
# opt$nodes = './whippet_nodes/Whippet_nodes.unambiguous_genes.gff'
# opt$breakpoints = 'isoSCM/compare_elav_pooled_all_wt_pooled_all.gtf'
# opt$min.confidence = 0.7
# opt$min.distance = 50
# opt$prefix = NULL

# params = list(min.confidence = 0.7, min.dist = 50, chrom.valid = c('2L','2R','3L','3R','4','X'))
params = list(min.confidence = opt$isoscm.confidence, 
              min.dist = opt$min.distance)

# utr3 as reference
txdb0 <- makeTxDbFromGFF(opt$gtf)
dm6.utr3 <- threeUTRsByTranscript(txdb0, use.names = TRUE)
dm6.tes <-  resize(transcripts(txdb0, use.names = TRUE), fix = 'end', width = 1)

dm6.ens96 <- import.gff(opt$gtf, feature.type = 'gene')
biotype.map <- deframe(mcols(dm6.ens96)[c('gene_id','gene_biotype')])

# Import and filter whippet nodes
breakpoints_all <- import.gff(opt$breakpoints)
whippet.nodes_all <- import.gff(opt$segments)

params$chrom.valid = as.character(unique(intersect(seqnames(whippet.nodes_all), seqnames(breakpoints_all))))

seqlevels(breakpoints_all, pruning.mode = 'coarse') <- params$chrom.valid
seqlevels(whippet.nodes_all, pruning.mode = 'coarse') <- params$chrom.valid


mcols(whippet.nodes_all)$is.utr3 <- overlapsAny(whippet.nodes_all, dm6.utr3)
whippet.nodes <- subset(whippet.nodes_all, 
                        biotype.map[Gene] %in% c('protein_coding','pseudogene') & 
                          seqnames %in% params$chrom.valid)


# Import and filter breakpoints
## Filter A) by: 
# - minimum confidence 
# - min distance
# - utr3 overlapping 
# - far enough from nearest TES
dd.tes <- distanceToNearest(breakpoints_all, dm6.tes)
mcols(breakpoints_all)$distance2tes <- Inf
mcols(breakpoints_all)$distance2tes[queryHits(dd.tes)] <- mcols(dd.tes)$distance

breakpoints <- subset(breakpoints_all, 
                        # as.double(confidence) > params$min.confidence & 
                        # type == 'changepoint' & 
                        overlapsAny(breakpoints_all, dm6.utr3) & 
                        distance2tes > params$min.dist)

if(length(breakpoints) == 0)
  stop('No breakpoints left after filtering.\nExit...')

# Filter B)
# require min distance from node start/end
breakpoints <- addDistanceInNode(breakpoints, whippet.nodes)

# d.end <- distanceToNearest(breakpoints, resize(whippet.nodes, fix = 'end', width = 1))
# d.start <- distanceToNearest(breakpoints, resize(whippet.nodes, fix = 'start', width = 1))
# 
# assertthat::assert_that(length(breakpoints) == length(d.start))
# assertthat::assert_that(all(queryHits(d.end) == queryHits(d.start)))
# 
# d.min <- apply(cbind(mcols(d.end)$distance, mcols(d.start)$distance), 1, min)
# mcols(breakpoints)$distance = d.min

breakpoints.filtered <- subset(breakpoints, distance > params$min.dist)

ii <- findOverlaps(breakpoints.filtered, whippet.nodes)
# Note: missed breakpoints are from ambiguous whippet nodes

# for each node, cut to pieces by breakpoints
cutNode <- function(start, end, node.template){
  node.template$width = NULL
    node.template$start = start
    node.template$end = end
  return(node.template)
}

cat('Breaking nodes\n')
whippet.mat <- as.data.frame(whippet.nodes)
whippet_split <- data.frame()
for(j in sort(unique(subjectHits(ii)))){
  # print(j)
  node.x <- whippet.mat[j,]
  i <- queryHits(subset(ii, subjectHits == j))
  breaks.start <- c(node.x$start, start(breakpoints.filtered[i]))
  breaks.end <- c(breaks.start[-1], node.x$end)
  
  node.y <- node.x
  nodes_new <- do.call(rbind, mapply(cutNode, breaks.start, breaks.end, MoreArgs = list(node.template=node.x), SIMPLIFY = FALSE))
  whippet_split <- rbind(whippet_split, nodes_new)
}

cat('Negative width nodes:\n')
whippet_split[which(whippet_split$end - whippet_split$start < 0),]

# whippet_split$width = whippet_split$end - whippet_split$start + 1

whippet_tmp <- makeGRangesFromDataFrame(whippet_split, keep.extra.columns = TRUE)
tmp = as.data.frame(whippet.nodes_all[!whippet.nodes_all %over% whippet_tmp])
tmp$width = NULL
whippet_novel <- sort(makeGRangesFromDataFrame(rbind(tmp,
                                       whippet_split), keep.extra.columns = TRUE))

tab0 <- as.data.frame(whippet_novel)
tab1 <- tab0 %>% 
  group_by(Gene) %>% 
  add_tally() %>%
  mutate(Node = ifelse(strand == '+', order(start), n - order(start) + 1))

tab.out <- data.frame(name = paste(tab1$Gene, tab1$Node, tab1$Type, sep = '_'), tab1)
assertthat::assert_that(all(table(tab.out$name) == 1 ), msg = 'Something went wrong. There are nodes in a gene with the same index. \nExit...')

# nodes to breaks
if(opt$debug){
  whippet2break <- whippet.nodes[unique(subjectHits(ii))]
  export.gff(whippet2break,paste0(opt$prefix,'.nodes_selected.gff'))
  # selected breakpoints
  breakpoints.selected <- breakpoints.filtered[queryHits(ii)]
  export.gff(breakpoints.selected,paste0(opt$prefix,'.breakpoints_selected.gff'))
}

write_tsv(tab.out, paste0(opt$prefix,'.saf'), col_names = FALSE)

whippet_novel <- makeGRangesFromDataFrame(tab.out, keep.extra.columns = TRUE)
export.gff3(whippet_novel, paste0(opt$prefix,'.gff'))

