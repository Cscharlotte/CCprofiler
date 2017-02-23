#' Filter a traces.object based on the Sibling Peptide Correlation (SibPepCorr).
#' @param traces.obj An object of type \code{traces.obj}.
#' @param protein_fdr_cutoff numeric specifying the FDR cutoff to be applied to select a sibling correlation cutoff
#' @param FFT numeric specifying the fraction of false targets (FFT). Default is 1.
#' @param absolute_spc_cutoff numeric specifying an absolute sibling correlation cutoff to be applied instead of FDR based cutoff estimation.
#' @param remove_decoys logical specifying if decoys should be removed. It is strongly recommended to keep the decoys. Default is FALSE.
#' @param plot logical TRUE or FALSE
#' @param PDF logical TRUE or FALSE
#' @param CSV logical TRUE or FALSE
#' @return traces.obj An object of type \code{traces.obj}.
#' @export

filterBySibPepCorr <- function(Traces,
                               protein_fdr_cutoff = 0.01,
                               FFT = 1,
                               absolute_spc_cutoff = NULL,
                               remove_decoys = TRUE,
                               plot = TRUE,
                               PDF = FALSE,
                               CSV = FALSE) {

  # If sibling peptide correlation hass not been calculated yet, do it
  if (!("SibPepCorr" %in% names(Traces$trace_annotation))){
    message("Sibling peptide correlation not yet calculated for this dataset\nCalculating SibPepCorr(spc)...")
    Traces <- calculateSibPepCorr(Traces, plot = plot, PDF = PDF)
  } else{
    message("Sibling peptide correlation values found...")
  }

  # get spc_cutoff from FDR estimation or direct
  decoys_contained = length(grep("^DECOY_", Traces$trace_annotation$protein_id)) > 0
  if (decoys_contained){
    message("Decoys found...\nEstimating FDR...")
    table <- ROC.SibPepCorr(Traces, FFT = FFT, plot = FALSE)
    cutoff_for_fdr <- roctable[FDR <= protein_fdr_cutoff, min(SibPepCorr_cutoff)]
    fdr_reached <- roctable[FDR <= protein_fdr_cutoff, max(FDR)]
    target_proteins_remaining <- roctable[FDR <= protein_fdr_cutoff, max(n_target_proteins)]
    message("Using SibCorrCutoff of ", round(cutoff_for_fdr, 4), "\nestimated FDR reached: ", round(fdr_reached, 4))
    spc_cutoff <- cutoff_for_fdr
  } else if (is.null(absolute_spc_cutoff)){
    message("No decoys found...\nPlease supply argument absolute_spc_cutoff\nor label decoy protein_ids with DECOY_ prefix")
  } else{
    message("Using absulte_spc_cutoff: ", absolute_spc_cutoff)
    spc_cutoff <- absolute_spc_cutoff
  }

  # subset Traces object
  Traces.filtered <- Traces
  Traces.filtered$trace_annotation <- Traces$trace_annotation[SibPepCorr >= spc_cutoff]
  Traces.filtered$traces <- Traces$traces[Traces$trace_annotation$SibPepCorr >= spc_cutoff]

  # compare target/decoy numbers before and after
  if (decoys_contained){
    # before
    Traces$trace_annotation$DECOY <- 0
    Traces$trace_annotation$DECOY[grep('^DECOY_', Traces$trace_annotation$protein_id)] <- 1
    target_proteins_before <- Traces$trace_annotation[DECOY == 0, length(unique(protein_id))]
    decoy_proteins_before <- Traces$trace_annotation[DECOY == 1, length(unique(protein_id))]
    protein_fdr_before = decoy_proteins_before/(decoy_proteins_before+target_proteins_before) * FFT
    # after
    Traces.filtered$trace_annotation$DECOY <- 0
    Traces.filtered$trace_annotation$DECOY[grep('^DECOY_', Traces.filtered$trace_annotation$protein_id)] <- 1
    target_proteins_after <- Traces.filtered$trace_annotation[DECOY == 0, length(unique(protein_id))]
    decoy_proteins_after <- Traces.filtered$trace_annotation[DECOY == 1, length(unique(protein_id))]
    protein_fdr_after = decoy_proteins_after/(decoy_proteins_after+target_proteins_after) * FFT
  } else{
    # before
    target_proteins_before <- Traces$trace_annotation[, length(unique(protein_id))]
    decoy_proteins_before <- 0
    protein_fdr_before <- NA
    # after
    target_proteins_after <- Traces.filtered$trace_annotation[, length(unique(protein_id))]
    decoy_proteins_after <- 0
    protein_fdr_after <- NA
  }

  # OUTPUT
  # Assemble small filter report table
  SibPepCorrFilter_report <- data.table(cutoff_placeholder = c("pre-filter", "post-filter", "fraction_removed:"),
             n_target_proteins = c(target_proteins_before, target_proteins_after, (1- target_proteins_after/target_proteins_before)),
             n_decoy_proteins = c(decoy_proteins_before, decoy_proteins_after, (1- decoy_proteins_after/decoy_proteins_before)),
             protein_fdr = c(protein_fdr_before, protein_fdr_after, ""))
  setnames(SibPepCorrFilter_report, "cutoff_placeholder", paste0("spc_cutoff: ", spc_cutoff))

  # CSV
  if (CSV){
    write.csv(roctable, file = "SibPepCorrFilter_IDFDRtable.csv", row.names = FALSE, quote = FALSE)
    write.csv(SibPepCorrFilter_report, file = "SibPepCorrFilter_report.csv", row.names = FALSE, quote = FALSE)
  }

  # PDF/ROC-like plot
  if (decoys_contained){
    if (PDF){
      pdf("SibPepCorrFilter_IDFDRplot.pdf")
    }
  if(plot){
	  par(mar=c(5.4,5.4,6.4,5.4))
      plot(roctable$FDR, roctable$n_target_proteins, type = "l", lty = 2,
           lwd = 2,
           main = paste0("ID-FDR Curve\n spc_cutoff: ", round(spc_cutoff, 4), "\ntarget proteins remaining: ",
                         target_proteins_remaining, "\nestimated FDR: ", round(fdr_reached, 4)),
           xlim = c(0,0.1), ylim = c(0, 1.02* max(roctable$n_target_proteins)), xlab='FDR',ylab='n')
      lines(roctable$FDR, roctable$n_true_target_proteins, lty = 1, lwd = 2)
      abline(h = target_proteins_remaining, lty = 2, col = "red")
      abline(v = fdr_reached, lty = 2, col = "red")
      points(fdr_reached, target_proteins_remaining, pty = 2, lwd = 3, col = "red")
      legend("bottomright", lty = c(2,1), legend = c("all target proteins", "true target proteins"))
    }
    if (PDF){
      dev.off()
    }
  }

  # remove decoys
  if (remove_decoys) {
    n_decoys <- length(grep("^DECOY", Traces.filtered$trace_annotation$protein_id))
    if ( n_decoys > 0){
      idx_decoys <- grep("^DECOY_",Traces.filtered$trace_annotation$protein_id)
      Traces.filtered$traces <- Traces.filtered$traces[-idx_decoys]
      Traces.filtered$trace_annotation<- Traces.filtered$trace_annotation[-idx_decoys]
      message(n_decoys, " decoys removed")
    } else {
      message("no decoys contained/removed")
    }
  }

  #return the filtered Traces data
  message("Proteins remaining in dataset: ", target_proteins_after)
  return(Traces.filtered)

}