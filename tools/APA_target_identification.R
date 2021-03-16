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
  make_option(c("-d", "--dexseq.table"), type="character",
              help = "DEX-seq table (tsv format)"),
  make_option(c("-g", "--gtf"), type="character",
              help="Reference annotation (ensembl gtf/gff)"),
  make_option(c("-s", "--segments"), type="character",
              help = "Segments, i.e. Whippet nodes (gtf/gff format)"),

  make_option(c("-w", "--width.min"), type="integer", default = 100,
              help = "Minimum segment size (nt) (default: %default)"),


  make_option("--padj.cutoff", type="double", default = 0.05,
              help = "DEXseq padj cutoff in node selection (default: %default)"),

  make_option("--write_locus", type="logical", action = 'store_true', default = FALSE,
              help = "Optional: Write locus file (gff)"),
  make_option("--DEXseq2gff", type="logical", action = 'store_true', default = FALSE,
              help = "Optional: Save DEX-seq as gff, includes coordinates (default: %default)"),
  make_option(c("--APA2gff"), type="logical", action = 'store_true', default = FALSE,
              help = "Output gff of output tsv)"),

  make_option(c("-o", "--outfileName"), type="character",
              help = "Name of file output file (tsv)")
)

opt <- parse_args(OptionParser(option_list=option_list))

print(opt)

suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tibble))


# get coordinates by id
nodes2granges <- function(dexseq.tab, whippet.nodes) {
  
  exon_id <- paste0(dexseq.tab$gene_id,'_',dexseq.tab$Node)
  mcols(whippet.nodes)$exon_id <- paste0(mcols(whippet.nodes)$Gene,'_',mcols(whippet.nodes)$Node)

  assertthat::assert_that(max(table(exon_id)) == 1, msg = 'Something went wrong with dexseq table ids')
  assertthat::assert_that(max(table(mcols(whippet.nodes)$exon_id)) == 1, msg = 'Something went wrong with whippet node ids')
  assertthat::assert_that(all(exon_id %in% mcols(whippet.nodes)$exon_id), msg = 'Not all dexseq ids are in whippet node ids')

  coords <- data.frame(row.names = mcols(whippet.nodes)$exon_id,
                       seqnames = seqnames(whippet.nodes),
                       start = start(whippet.nodes),
                       end = end(whippet.nodes),strand = strand(whippet.nodes))
  gr <- makeGRangesFromDataFrame(data.frame(coords[exon_id,], dexseq.tab), keep.extra.columns = TRUE)
  
  return(gr)
}

isUTR3_and_sameGene <- function(exonset, txdb, tx2gid.map){
  
  is.utr3_and_sameGene = rep(FALSE, length(exonset))
  
  utr3.df <- as.data.frame(threeUTRsByTranscript(txdb, use.names = TRUE))
  utr3.df$gene_id <- tx2gid.map[utr3.df$group_name]
  
  utr3.gr <- makeGRangesFromDataFrame(utr3.df, keep.extra.columns = TRUE)
  
  
  q.overlaps <- findOverlaps(query = exonset, subject = utr3.gr)
  gid.dexseq = mcols(exonset)$gene_id[queryHits(q.overlaps)]
  gid.utr3 = mcols(utr3.gr)$gene_id[subjectHits(q.overlaps)]
  is.utr3_and_sameGene[queryHits(q.overlaps)] <-  gid.dexseq == gid.utr3
  
  return(is.utr3_and_sameGene)
}

exon.crawler <- function(exonset, padj.cutoff){
  
  cols.required = c('Node','is.anchor','is.utr3','log2fold_cond_ctrl','padj')
  assertthat::assert_that(all(cols.required %in% colnames(exonset)))
  
  if(nrow(exonset) < 2 | 
     !any(exonset$padj < padj.cutoff & exonset$is.utr3, na.rm = TRUE) | 
     sum(exonset$is.utr3,na.rm = TRUE) < 2)
    return(NULL)
  
  exonset <- exonset[order(exonset$Node, decreasing = TRUE, na.last = TRUE),]
  exonset$uid <- 1:nrow(exonset)
  exonset.utr3 <- subset(exonset, is.utr3)
  
  j <- which(exonset.utr3$is.anchor)
  if(length(j) != 1 | all(is.na(exonset.utr3$padj))){
    cat('Omit exonset.utr3\n')
    print(exonset.utr3 %>% as.data.frame)
    return(NULL)
  }
  
  lfc.direction <- sign(exonset.utr3$log2fold_cond_ctrl[j])
  # strategy 1) check padj(i+1)
  
  # startegy 2) two-pass approach
  i = j 
  while(i <= nrow(exonset.utr3) &
        !(is.na(exonset.utr3$padj[i])) &
        exonset.utr3$padj[i] < padj.cutoff &
        sign(exonset.utr3$log2fold_cond_ctrl[i]) == lfc.direction &
        exonset.utr3$is.utr3[i]){
    i = i + 1
  }
  
  
  if(i > nrow(exonset.utr3)){
    if(j == i-1)
      return(exonset[exonset.utr3$uid[j],,drop=FALSE] %>% mutate(is.firstUTR3 = TRUE))
    
    return(exonset[exonset.utr3$uid[(i-1):j],,drop=FALSE] %>% mutate(is.firstUTR3 = TRUE))
  } else {
    return(exonset[exonset.utr3$uid[i:j],,drop=FALSE] %>% mutate(is.firstUTR3 = FALSE))
  }
}

run_exon.crawler = function(exonset, padj.cutoff){
  
  t.list = split(exonset, f = exonset$gene_id)
  t.list = lapply(t.list, function(x) x[order(x$Node, decreasing = TRUE, na.last = TRUE),])
  
  nodes.diff_exons_list <- lapply(t.list, exon.crawler, padj.cutoff = padj.cutoff)
  
  nodes.apa <- do.call(rbind, nodes.diff_exons_list)
  return(nodes.apa)
}

aggregateDEXseqWeightedMean <- function(apa_nodeset, cat.cols, metric.cols) {
  
  apa_nodeset = as.data.frame(apa_nodeset)
  
  paste1=function(x) paste0(unique(x), collapse ='_')
  
  i <- which.min(apa_nodeset$Node)
  breakpoint = apa_nodeset[i,]
  
  if(breakpoint$is.firstUTR3) {
    if(breakpoint$strand == '+')   
      breakpoint$end = breakpoint$start 
    if(breakpoint$strand == '-')
      breakpoint$start = breakpoint$end 
  } else {
    if(breakpoint$strand == '+')
      breakpoint$start = breakpoint$end 
    if(breakpoint$strand == '-') 
      breakpoint$end = breakpoint$start
    
    if(nrow(apa_nodeset)-1 < 1)
      print(apa_nodeset)
    
    apa_nodeset = apa_nodeset[-i,]
    assertthat::assert_that(nrow(apa_nodeset) > 0, msg = 'You found an empty apa_nodeset in \'not firstUTR3\'')
  }
  
  padj.i <- which.min(apa_nodeset$padj)
  for(n0 in metric.cols) breakpoint[,n0] <- apa_nodeset[padj.i,n0]
  for(n0 in cat.cols) breakpoint[,n0] <- paste1(apa_nodeset[,n0])
  
  #change log2FC to weighted average of all nodes
  breakpoint$log2fold_cond_ctrl = weighted.mean(apa_nodeset$log2fold_cond_ctrl, 
                                                apa_nodeset$width/sum(apa_nodeset$width), na.rm = TRUE)
  # mean of normalized exonBaseMean
  breakpoint$exonBaseMean = mean(apa_nodeset$exonBaseMean/apa_nodeset$width, na.rm = TRUE) * 1e3
  
  breakpoint$node_most_significant <- apa_nodeset$Node[padj.i]
  breakpoint$n.nodes = nrow(apa_nodeset)
  breakpoint$width = sum(apa_nodeset$width)
  
  col_select = c('seqnames','start','end','width','strand', 
                 c(cat.cols, metric.cols), 
                 'node_most_significant','n.nodes', 'is.firstUTR3')
  return(breakpoint[col_select])
}

mergeLocusAndBreakpoint <- function(breakpoints, exonset){
  
  .makeSet = function(breakpoint, exonset){
    
    nodeset <- as.integer(unlist(str_split(breakpoint$Node,'_')))
    
    targets <- subset(exonset, gene_id %in% breakpoint$gene_id & Node %in% nodeset)
    start.locus = min(targets$start, na.rm = TRUE)
    end.locus = max(targets$end, na.rm = TRUE) 
    
    breakpoint$type = 'breakpoint'
    breakpoint$source = 'APA_targets_identification'
    
    locus <- data.frame(seqnames = unique(targets$seqnames),
                        start = c(breakpoint$start, start.locus),
                        end = c(breakpoint$end, end.locus),
                        strand = breakpoint$strand,
                        type = c('breakpoint','locus'),
                        source = 'APA_targets_identification',
                        gene_id = unique(targets$gene_id),
                        gene_name = unique(targets$gene_name),
                        stringsAsFactors = FALSE)
    return(locus)
  }
  
  breakpoint.set = split(breakpoints, breakpoints$gene_id)
  node_set = do.call(rbind, lapply(breakpoint.set, .makeSet, exonset))
  
  return(node_set)
}

# load data. txdb (for 3'utrs), DEXseq table, whippet nodes (for coordinates)
anno.ref <- import.gff(opt$gtf, feature.type = 'gene')
anno.transcripts <- import.gff(opt$gtf, feature.type = 'transcript')
tx2gid = deframe(mcols(anno.transcripts)[c('transcript_id','gene_id')])

gid2symbol <- deframe(mcols(anno.ref)[c('gene_id','gene_name')])

txdb0 <- makeTxDbFromGFF(opt$gtf)
anno.utr3 <- threeUTRsByTranscript(txdb0, use.names = TRUE)

whippet.nodes <- import.gff(opt$segments)

dexseq_tab = read_tsv(opt$dexseq.table)

dexseq_tab$Node <- as.integer(gsub('^E','',dexseq_tab$featureID))
colnames(dexseq_tab)[grepl('groupID',colnames(dexseq_tab))] = 'gene_id'
dexseq_tab$gene_name <- gid2symbol[dexseq_tab$gene_id]

# add coordinates ,filter for UTR3
gr0 <- nodes2granges(dexseq_tab, whippet.nodes)
mcols(gr0)$is.utr3 = isUTR3_and_sameGene(gr0, txdb0, tx2gid.map = tx2gid)
dexseq_tab <- as.data.frame(gr0)


params = list()
params$padj.cutoff = opt$padj.cutoff

# Identify last significant nodes
nodes.sign = subset(dexseq_tab, padj < params$padj.cutoff & is.utr3)

nodes.sign_last = nodes.sign %>% 
  dplyr::group_by(gene_id) %>% 
  dplyr::filter(!is.na(log2fold_cond_ctrl)) %>%
  dplyr::slice(which.max(Node))

n = paste0(nodes.sign_last$gene_id,'_',nodes.sign_last$Node )
apa.genes <- unique(nodes.sign_last$gene_id);length(apa.genes)

# add annotation to original dexseq table
dexseq_annotated = dexseq_tab %>% 
  dplyr::group_by(gene_id) %>% 
  dplyr::mutate(node.id = paste0(gene_id,'_',Node)) %>% 
  dplyr::mutate(is.anchor = node.id %in% n)

dexseq_apa.genes = dexseq_annotated %>%
  filter(gene_id %in% apa.genes)

dexseq_apa.genes_filtered = dexseq_apa.genes %>% 
  filter(width > opt$width.min) %>% 
  group_by(gene_id) %>% 
  filter(any(is.anchor))

nodes.apa = run_exon.crawler(dexseq_apa.genes_filtered, params$padj.cutoff)

# aggregate values across changed exons

nodes.apa_list <- split(nodes.apa, nodes.apa$gene_id)
res <- do.call(rbind, 
               lapply(nodes.apa_list,
                 aggregateDEXseqWeightedMean,
                 cat.cols = c('Node','gene_id','gene_name'), 
                 metric.cols = c('exonBaseMean','dispersion','stat','pvalue','padj','ctrl','cond','log2fold_cond_ctrl')
               ))

write_tsv(res, opt$outfileName)

if(opt$DEXseq2gff){
  dexseq.gff_filename = paste0(dirname(opt$outfileName),'/', gsub('.tsv','.gff',basename(opt$dexseq.table)))
  gr1 <- makeGRangesFromDataFrame(dexseq_annotated %>% 
                                    select(-contains('_pooled_'),-contains('is.utr3'),-contains('is.anchor'),-contains('node.id')), 
                                           keep.extra.columns = TRUE)
  export.gff3(gr1, dexseq.gff_filename)
}

if(opt$APA2gff){
  apa.gff_filename = paste0(gsub('\\.tsv$','',opt$outfileName),'.gff')

  gr0 <- makeGRangesFromDataFrame(res, keep.extra.columns = TRUE)
  export.gff3(gr0, apa.gff_filename)
}

if(opt$write_locus){
  locus_filename =  paste0(gsub('\\.tsv$','',opt$outfileName),'.locus.gff')
  
  locusset = mergeLocusAndBreakpoint(res, dexseq_annotated)
  
  nodeSet0 <- makeGRangesFromDataFrame(locusset, keep.extra.columns = TRUE)
  export.gff3(nodeSet0, locus_filename)
}
cat('Done\n')
