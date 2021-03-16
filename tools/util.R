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
