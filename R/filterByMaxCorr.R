
getCorrelationMatrices <- function(traces){
  genes <- unique(traces$trace_annotation$protein_id)
  intMat <- getIntensityMatrix(traces)
  genePepList <- lapply(genes, FUN = function(gene){
    peps <- traces$trace_annotation[protein_id == gene, id]
    res <- intMat[peps,]
    return(t(res))
  })
  names(genePepList) <- genes
### find maximum corr for every peptide in each protein
  corrMatrices <- lapply(genePepList, function(gene){
    genecorr <- cor(gene)
    # genecorr <- cor(gene, method = "spearman")
    return(genecorr)
  })
  names(corrMatrices) <- names(genePepList)
  return(corrMatrices)
}

filterByMaxCorr <- function(traces, corrMatrices = NULL, cutoff = 0.85, plot = FALSE) {
  if(is.null(corrMatrices)){
    corrMatrices <- getCorrelationMatrices(traces)
  }
  maxCorrMatrices <- lapply(corrMatrices, function(gene){
    genecorr <- gene
    genecorr[genecorr == 1] <- NA
    maxcorr <- rowMaxs(genecorr, na.rm = T)
    names(maxcorr) <- rownames(genecorr)
    return(maxcorr)
  })
  
  names(maxCorrMatrices) <- names(genePepList)
  allMaxCorrs <- unlist(maxCorrMatrices)
  names(allMaxCorrs) <- gsub(".*\\.", "", names(allMaxCorrs))
  if (plot == TRUE){
    hist(allMaxCorrs, breaks = 30)
    plot(density(allMaxCorrs), xlim = c(-0.2, 0.5))
  }
  
  filterpeps <- names(allMaxCorrs)[allMaxCorrs > cutoff]
  
  traces_filt <- subset(traces, filterpeps, "id")
  return(traces_filt)
}