suppressPackageStartupMessages(library("optparse"))

# specify our desired options in a list
# by default OptionParser will add an help option equivalent to 
# make_option(c("-h", "--help"), action="store_true", default=FALSE, 
#               help="Show this help message and exit")
option_list <- list( 
  # make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
  #             help="Print extra output [default]"),
  # make_option(c("-q", "--quietly"), action="store_false", 
  #             dest="verbose", help="Print little output"),
  make_option(c("-t", "--table"), type="character", 
              help="featureCount table. Requires rownames to be from Whippet node annotation <Gene>_<Node>_<Type>"),
  make_option(c("-s", "--samplesheet"), type="character",
              help = "Sample sheet to passed to DEXseq. Needs to match 'basename' of column names of -t"),
  make_option(c("-o", "--outfileName"), type="character",
              help = "Name for output table (tsv file)")
)

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
opt <- parse_args(OptionParser(option_list=option_list))
print(opt)

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(DEXSeq))
suppressPackageStartupMessages(library(janitor))

runDEXseq <- function(node.tab, sampleData) {
  # sampleData <- data.frame(sample = colnames(NodeCountsElav_TE[,c(2:5)]),
  #                          condition=gsub("(.*)_[12]","\\1",colnames(NodeCountsElav_TE[,c(2:5)])))
  require(DEXSeq)
  
  groupID = gsub(":E.*","",rownames(node.tab))
  featureID = gsub(".*:(E.*)","\\1",rownames(node.tab))
  
  assertthat::assert_that(all(colnames(sampleData) %in% c('name','condition')))
  
  dxd <- DEXSeqDataSet( node.tab[,sampleData$name], sampleData, 
                        design= ~ name + exon + condition:exon, 
                        featureID, groupID, featureRanges=NULL, 
                        transcripts=NULL, alternativeCountData=NULL)
  BPPARAM = MulticoreParam(workers=4)
  dxd = estimateSizeFactors( dxd )
  dxd = estimateDispersions( dxd, BPPARAM=BPPARAM)
  # plotDispEsts( dxd )
  
  dxd = testForDEU( dxd, BPPARAM=BPPARAM)
  dxd = estimateExonFoldChanges(dxd, BPPARAM=BPPARAM)
  
  dxr1 = DEXSeqResults( dxd )
  # plotMA( dxr1, cex=0.8 )
  
  sampleAnnotation(dxd)
  
  return(dxr1)
}


# whippet.featurects <- read_tsv('./Whippet_nodes.featureCounts.tsv', comment ='#')
# whippet.featurects <- read_tsv('./Whippet_nodes.split.featureCounts.tsv', comment ='#')

whippet.featurects <- read_tsv(opt$table, comment = '#')
colnames(whippet.featurects) <- gsub('.bam','',basename(colnames(whippet.featurects)))
whippet.featurects <- as.data.frame(clean_names(whippet.featurects))

samplesheet = as.data.frame(read_tsv(opt$samplesheet))
samplesheet$name <- make_clean_names(samplesheet$name)
samplesheet$condition <- factor(samplesheet$condition, levels = unique(samplesheet$condition))

assertthat::assert_that(all(samplesheet$name %in% colnames(whippet.featurects)), msg = 'Mismatching names in samplesheet')

# t0 <- grep(pattern,colnames(whippet.featurects), value=TRUE)
# sample.table <- data.frame(sample = t0, condition = as.factor(ifelse(grepl('wt',t0),'control','condition')))
# sample.table$condition <- relevel(sample.table$condition, 'control')

t0 <- str_split_fixed(whippet.featurects$geneid, '_', n = 3)
rownames(whippet.featurects) <- paste0(t0[,1],':E',t0[,2])

# parse data set
node.tab <- whippet.featurects[,samplesheet$name]

# run DEXseq
dexseq0 <- runDEXseq(node.tab = round(node.tab), sampleData = samplesheet)

# write_tsv(as.data.frame(dexseq0), './results/DEXseq.wt_elav.whitelist.alpha.tsv')
# write_tsv(as.data.frame(dexseq0), './results/DEXseq.wt_elav.whippet_split.whitelist.alpha.tsv')

tab0 = as.data.frame(dexseq0)
tab0$genomicData = NULL
colnames(tab0)[grepl('^countData',colnames(tab0))] = samplesheet$name

write_tsv(tab0, opt$outfileName)

