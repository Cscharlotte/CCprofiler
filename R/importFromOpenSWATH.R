# Due to: http://stackoverflow.com/questions/24501245/data-table-throws-object-not-found-error
.datatable.aware=TRUE

#' Import peptide profiles from an OpenSWATH experiment.
#' @description This is a convenience function to directly import peptide profiles
#'  from an OpenSWATH experiment (after TRIC alignment). The peptide intensities are calculated by summing all charge
#'  states. Alternativley the MS1 signal can be used for quantification.
#' @import data.table
#' @param data Quantitative MS data in form of OpenSWATH result file or R data.table.
#' @param annotation_table file or data.table containing columns `filename` and
#'     `fraction_number` that map the file names (occuring in input table filename column)
#'     `to a SEC elution fraction.
#' @param rm_requantified Logical, whether requantified (noise) peak group quantities
#'     (as indicated by m_score = 2) should be removed, defaults to \code{TRUE}.
#' @param rm_decoys Logical, whether decoys should be removed, defaults to \code{FALSE}.
#' @param MS1Quant Logical, whether MS1 quantification should be used, defaults to \code{FALSE}.
#' @param verbose Logical, whether to print progress message into console, defaults to \code{TRUE}.
#' @return An object of class Traces containing
#'     "traces", "traces_type", "traces_annotation" and "fraction_annotation" list entries
#'     that can be processed with the herein contained functions.
#' @export
#' @examples
#'   input_data <- exampleOpenSWATHinput
#'   annotation <- exampleFractionAnnotation
#'   traces <- importFromOpenSWATH(data = input_data,
#'                                 annotation_table = annotation,
#'                                 rm_requantified = TRUE)
#'   summary(traces)
#'

importFromOpenSWATH <- function(data,
                                annotation_table,
                                rm_requantified=TRUE,
                                rm_decoys = FALSE,
                                MS1Quant=FALSE,
                                verbose=TRUE){

  ## test arguments
  if (missing(data)){
        stop("Need to specify data in form of OpenSWATH result file or R data.table.")
  }
  if (missing(annotation_table)){
        stop("Need to specify annotation_table.")
  }

  ## read data & annotation table

  if (class(data)[1] == "character") {
    if (file.exists(data)) {
      message('reading results file ...')
      data  <- data.table::fread(data, header=TRUE)
    } else {
      stop("data file doesn't exist")
    }
  } else if (all(class(data) != c("data.table","data.frame"))) {
    stop("data input is neither file name or data.table")
  }

  if (class(annotation_table)[1] == "character") {
    if (file.exists(annotation_table)) {
      message('reading annotation table ...')
      annotation_table  <- data.table::fread(annotation_table)
    } else {
      stop("annotation_table file doesn't exist")
    }
  } else if (all(class(annotation_table) != c("data.table","data.frame"))) {
    stop("annotation_table input is neither file name or data.table")
  }

  ## remove non-proteotypic discarding/keeping Decoys
  message('removing non-unique proteins only keeping proteotypic peptides ...')
  if (rm_decoys == TRUE){
    message('remove decoys ...')
    data <- data[grep("^1/", data$ProteinName)]
  } else {
    data <- data[c(grep("^1/", data$ProteinName), grep("^DECOY_1/", data$ProteinName))]
  }
  data$ProteinName <- gsub("1\\/","",data$ProteinName)


  ## convert ProteinName to uniprot ids
  if (length(grep("\\|",data$ProteinName)) > 0) {
    message('converting ProteinName to uniprot ids ...')
    decoy_idx <- grep("^DECOY_",data$ProteinName)
    data$ProteinName <- gsub(".*\\|(.*?)\\|.*", "\\1", data$ProteinName)
    data$ProteinName[decoy_idx] <- paste0("DECOY_",data$ProteinName[decoy_idx])
    # the above does not work further downstream because "1/" is removed
    #data$ProteinName <- extractIdsFromFastaHeader(data$ProteinName)
  }

  ## subset data to some important columns to save RAM
  if (MS1Quant == TRUE) {
  column_names <- c('transition_group_id', 'ProteinName','FullPeptideName',
                    'filename', 'Sequence', 'decoy', 'aggr_prec_Peak_Area',
                    'd_score', 'm_score')
  }else{
  column_names <- c('transition_group_id', 'ProteinName','FullPeptideName',
                      'filename', 'Sequence', 'decoy', 'd_score', 'm_score', 'Intensity')
  }
  data_s <- subset(data, select=column_names)

  ## Use aggregated precursor (MS1) area as Intensity column if selected
  if (MS1Quant == TRUE){
    setnames(data_s, 'aggr_prec_Peak_Area', 'Intensity')
  }

  rm(data)
  gc()

  ## remove requantified values if selected
  if (rm_requantified == TRUE) {
      data_s <- data_s[m_score < 2, ]
  }

  ## add fraction number column to main dataframe
  fraction_number <- integer(length=nrow(data_s))
  files <- annotation_table$filename
  data_filenames <- data_s$filename

  if (length(files) != length(unique(data_filenames))) {
      stop("Number of file names in annotation_table does not match data")
  }


  for (i in seq_along(files)) {
      idxs <- grep(files[i], data_filenames)
      fraction_number[idxs] <- annotation_table$fraction_number[i]
      if (verbose) {
        message(paste("PROCESSED", i, "/", length(files), "filenames"))
      }
  }
  data_s <- cbind(data_s, fraction_number)

  ## Assemble and output result "traces" object
  traces_wide <-
      data.table(dcast(data_s, ProteinName + FullPeptideName ~ fraction_number,
                       value.var="Intensity",
                       fun.aggregate=sum))

  traces_annotation <- data.table(traces_wide[,c("FullPeptideName", "ProteinName"), with = FALSE])
  setcolorder(traces_annotation, c("FullPeptideName", "ProteinName"))
  setnames(traces_annotation,c("FullPeptideName", "ProteinName"),c("id","protein_id"))

  traces <- subset(traces_wide, select = -ProteinName)
  traces[,id:=FullPeptideName]
  traces[,FullPeptideName:=NULL]

  nfractions <- ncol(traces)-1
  fractions <- as.numeric(c(1:nfractions))
  fraction_annotation <- data.table(id=fractions)

  traces_type = "peptide"

  result <- list("traces" = traces,
                 "trace_type" = traces_type,
                 "trace_annotation" = traces_annotation,
                 "fraction_annotation" = fraction_annotation)
  class(result) <- "traces"
  names(result) <- c("traces", "trace_type", "trace_annotation", "fraction_annotation")

  # sort both trace_annotation and traces for consistency (also done in annotateTraces function)
  setorder(result$trace_annotation, id)
  setorder(result$traces, id)
  .tracesTest(result)
  return(result)
}
