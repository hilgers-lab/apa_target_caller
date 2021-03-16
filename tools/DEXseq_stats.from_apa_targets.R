# Given a set of target segments, extract DEXseq measurements according to segmentation. 
# Note, segmentation must be identical to target segments
#
# Speciality in this script: Segmentations are from various sources. 
# 
library(optparse)

# segments.gff 
# apa_targets.tsv
# dexseq.tsv # must match segments.gff

option_list <- list( 
  make_option(c("-s", "--segmentation"), type='character',
              help="Segments used for DEXSeq quantification [gff] (Segments_split.gff of UTR3_quantification)"),
  make_option(c("-d", "--dexseq"), type='character',
              help="dexseq table. Quantification based on --segments [tsv]"),
  make_option(c("-t", "--target_coordinates"), type='character',
              help="APA segments to extract DEXseq measures for [tsv]"),
  
  make_option(c("-o", "--outprefix"), type='character',
              help="Output table, summary stats per gene [tsv]"),
  
  
  make_option("--per_node", action='store_true',default=FALSE,
              help="create output also per node"),

  make_option("--gff", action='store_true',default=FALSE,
              help="create gff from output table(s)")
)


# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
opt <- parse_args(OptionParser(option_list=option_list))

if(FALSE){
  opt$segmentation = "results_early/nfne_wt.early/Annotation/Segments_split.gff"
  opt$dexseq="results_early/nfne_wt.early/DEXseq/Project_1384.nFne.early.segments_split.dexseq.tsv"
  opt$target_coordinates="Table_S1_APA_target_genes.elav_early.segments.tsv"
  opt$outprefix="test_run"
  opt$per_node=TRUE
  opt$gff=TRUE
}

segmentation_f = ifelse(file.exists(opt$segmentation), 
                      opt$segmentation,
                      stop('segmentation file not found'))
dexseq_f = ifelse(file.exists(opt$dexseq),
                opt$dexseq, 
                stop('dexseq table does not exist'))
# apa_targets required columns: Nodes, gene_id
apa_targets_f = ifelse(file.exists(opt$target_coordinates),
                     opt$target_coordinates,
                     stop('APA target file not found'))
                     

suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(rtracklayer))
# suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tibble))

expand.data.frame = function(df, group, sep){
  group.split = strsplit(group, sep)
  i = rep(1:nrow(df), times = elementNROWS(group.split))
  dd = data.frame(df[i,], group = as.character(unlist(group.split)))
  return(dd)
}

cat(">>> Importing Tables\n")

segmentations = import.gff(segmentation_f) %>% as.data.frame()

apa_targets = read_tsv(apa_targets_f)
apa_targets = apa_targets %>% mutate(Nodeset = as.character(apa_targets$Nodeset))

dexseq.tab = read_tsv(dexseq_f)
dexseq.tab$Node = gsub('^E','',dexseq.tab$featureID)
dexseq.tab$gene_id = dexseq.tab$groupID

# expand apa_targets.per_node by segment

cat(">>> Expanding targets per node\n")
apa_targets.per_node = expand.data.frame(apa_targets, apa_targets$Nodeset, '_')
apa_targets.per_node = apa_targets.per_node %>% 
  dplyr::rename(Node = group)

## Add node width annotation from segmentation
# tmp_set = list()
# for(type in c('early','late','neuronal')){
#   tmp_set[[type]] = apa_targets.per_node %>% 
#     filter(segmentation == type) %>% 
#     left_join(segmentations[[type]] %>% dplyr::select(Gene, Node, width), 
#               by = c('gene_id' = 'Gene','Node')) %>% 
#     mutate(width = width.y) %>% 
#     dplyr::select(-width.x,-width.y)
# }
# apa_targets.per_node = do.call(rbind, tmp_set)

cat(">>> Collecting annotation targets per node\n")
tmp = apa_targets.per_node%>% dplyr::select(-width) %>% 
  left_join(segmentations %>% dplyr::select(Gene, Node, width), by = c('gene_id' = 'Gene', 'Node')); tmp %>% dim
apa_targets.per_node = tmp
## Add DEXseq measures by gene_id, Node
# targets = list()
# for(type in c('early','late','neuronal')){
#   targets[[type]] = apa_targets.per_node %>%
#     filter(segmentation == type) %>% 
#     left_join(dexseq_mutant.set[[type]] %>%
#                 dplyr::select(gene_id, Node, exonBaseMean, log2fold_cond_ctrl, padj), 
#               by = c('gene_id','Node'))
# }
# tab.per_node = do.call(rbind,targets)
cat(">> Adding DEXseq measures\n")
tab.per_node = apa_targets.per_node %>%
  left_join(dexseq.tab %>%
              dplyr::select(gene_id, Node, exonBaseMean, log2fold_cond_ctrl, padj), 
            by = c('gene_id','Node'))

cat(">>> Summarize measure per gene\n")
# Summarize per gene_id
tab.per_gene = tab.per_node %>% 
  group_by(gene_id, seqnames, start, end) %>% 
  summarize(nodeset = paste(unique(Node), collapse = '_'),
            # min padj
            padj.min = ifelse(all(is.na(padj)), NA, min(padj, na.rm = TRUE)),
            # average exonBaseMean
            exonBaseMean.mean = mean(exonBaseMean, na.rm = TRUE), 
            # log2fold_cond_ctrl
            log2fold_cond_ctrl.weighted = weighted.mean(log2fold_cond_ctrl, width/sum(width), na.rm = TRUE),
            # weighed mean of log2FoldChange, weighted by width proportion
            exonBaseMean.weighted = mean(exonBaseMean/width, na.rm = TRUE) * 1e3,
            node_index_most_significant = ifelse(all(is.na(padj)), NA, which.min(padj)))


outprefix = opt$outprefix
cat(">>> Writing files\n")
write_tsv(tab.per_gene,  paste0(outprefix,'.apa_targets.gene.tsv'))

if(opt$gff) {
  library(GenomicFeatures)
}

if(opt$gff){
  gff0 = tab.per_gene %>% makeGRangesFromDataFrame(keep.extra.columns = TRUE)
  names(gff0) = paste(mcols(gff0)$gene_id, mcols(gff0)$nodeset, sep = '.')
  export.gff2(gff0,  paste0(outprefix,'.apa_targets.gene.gff'))
}

if(opt$per_node) { # requires mapping of coordinates to each node, instead of whole segment. respect hierarchy
  # tmp_set = list()
  # for(type in c('early','late','neuronal')){
  #   tmp_set[[type]] = tab.per_node %>% dplyr::select(-seqnames,-start,-end,-strand) %>%
  #     filter(segmentation == type) %>%
  #     left_join(segmentations[[type]] %>% dplyr::select(Gene,'seqnames','start','end', 'strand', 'Node'),
  #               by = c('gene_id' = 'Gene','Node'))
  # }
  # tab.per_node.coords = do.call(rbind, tmp_set)
  
  tab.per_node.coords = tab.per_node %>% dplyr::select(-seqnames,-start,-end,-strand) %>%
        left_join(segmentations %>% dplyr::select(Gene,'seqnames','start','end', 'strand', 'Node'),
                  by = c('gene_id' = 'Gene','Node'))
  write_tsv(tab.per_node.coords, paste0(outprefix,'.apa_targets.nodes.tsv'))
  
  if(opt$gff) {
    gff0 = tab.per_node.coords %>% makeGRangesFromDataFrame(keep.extra.columns = TRUE)
    names(gff0) = paste(mcols(gff0)$gene_id, mcols(gff0)$Node,sep = '_')
    export.gff2(gff0, paste0(outprefix,'.apa_targets.nodes.gff'))
  }
}
