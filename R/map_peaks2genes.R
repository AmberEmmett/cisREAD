#' Map ATAC-seq peaks to nearby genes

#' This function maps accessible genomic regions to protein coding genes 
#' within a given distance.

#' @param peaks An ATAC-seq count matrix, where the first 6 columns are in BED 
#' format, and columns afterwards are samples. Each row of the matrix should 
#' give normalised, log2 transformed ATAC-seq counts for each ATAC-seq peak.
#' @param genome The genome annotation, only "hg38" is currently supported 
#' @param distance The maximum genomic distance, in base pairs, allowed between a 
#' gene's TSS and the centre of the ATAC-seq peak (default = 100000)
#' @return A dataframe where each peak is mapped to all genes within
#' the given distance
#' @example 
#' head(ATAC_peaks)
#' peaks2genes <- map_peaks2genes(peaks = ATAC_peaks, genome = "hg38", distance = 100000)
#' head(peaks2genes)
#' @export
map_peaks2genes <- function(peaks, genome = "hg38", distance=100000) {

  genome <- rlang::arg_match(genome)

  if (genome == "hg38") {
    genes = hg38
  }

  genes$upper <- genes$TSS + distance
  genes$lower <- genes$TSS- distance
  
  names(peaks)[1:6] <- c("Chr", "Start", "End", "Peak", "width", "strand")
  nCells <- length(peaks) - 6

  peaks$middle <- round((peaks$Start + peaks$End) / 2)
  peaks2genes <-  sqldf::sqldf("select a.*, b.*
      from peaks a left join genes b
      on a.Chr = b.chr and a.middle between
              b.lower and
              b.upper")

  peaks2genes <- na.omit(peaks2genes)
  peaks2genes$Distance <- peaks2genes$middle - peaks2genes$TSS
  peaks2genes$Distance <- ifelse(peaks2genes$Strand == 1, peaks2genes$Distance, (peaks2genes$Distance * -1))
  peaks2genes$Strand <- ifelse(peaks2genes$Strand == 1, "+", "-")
  peaks2genes <- peaks2genes[ ,c(1:4, 
                                 (nCells + 8):(nCells + 9),
                                 (nCells + 11):(nCells + 12),
                                 (nCells + 15), 7:(nCells+6))]
  peaks2genes
}
