#' Performs the cisREAD workflow in one step upon a list of genes, from ATAC-seq , RNA-seq and TF footprint inputs

#' This function is a wrapper of map_peaks2genes, find_coCREs and select_coCREs which can be iterated over a list
#' of genes.

##' @param peaks An ATAC-seq count matrix, where the first 6 columns are in BED 
#' format, and columns afterwards are samples. Each row of the matrix should 
#' give normalised, log2 transformed ATAC-seq counts for each ATAC-seq peak.
#' @param TFprofile A binary matrix where each row is an ATAC-seq peak, and each column indicates the presence (1)
#' or absence (0) of a TF footprint in a given differentiation stage
#' @param RNA A gene expression matrix where each row is a gene's ENSEMBL ID (if not then please convert first)
#' and each column gives the gene's expression (log2 normalised counts) in a given sample
#' @param genes A character vector of protein-coding gene symbols to find and select coCREs for - 
#' e.g. c("PRDM1", "MYC", "IRF4","BCL6")
#' @param genome The genome annotation, only "hg38" is currently supported 
#' @param distance The maximum genomic distance allowed between a candidate CRE and the gene's TSS. Cannot be 
#' larger than the distance used in maps_peaks2genes
#' @param minTFevents Used to filter out candidate CREs not bound by any TFs. All candidate CREs must have at least
#' this number of 1's in the TFprofile matrix (default 2)
#' @param coCREcutoff The minimum similarity score (TF footprint dice similarity * ATAC-seq pearson correlation)
#' used to draw an edge between two candidate CREs in network construction. Should be a value between 0 and 1
#' (default 0.3), ower values will result in looser communities and higher values in tighter communities
#' @param weights A vector giving the weights assigned to TF similarity (1st value) and chromatin accessibility
#' correlation (2nd value). The defualt - c(1,1) - gives equal weighting to both, using c(0,1) would perform
#' community detection solely on chromatin accessibility
#' @param lambda The value of lambda used for model selection during cross-validation, one of 
#' "lambda_min" (default) and "lambda_se". "lambda_min" will construct the model with the mean
#' minimum cross-validated squared error, "lambda_se" will construct the model within 1 standard error 
#' of "lambda_min", and will result in a sparser model
#' @param pval The p-value alpha level for assigning statistical significance. coCRE predictors with
#' p < this value will be considered statistically significant
#' @param padj The Benjamini-Hockberg adjusted p-value alpha level for assigning statistical significance to
#' LASSO models. Allows for multiple testing correction when applying cisREAD to many genes.
#' (co)CRE predictors with p < pval, from gene-specific models with p < padj will be considered statistically
#' significant
#' @return A cisREADResult object, listing significant, selected and all candidate CRE communities
#' @example 
#' #Example input data
#' head(ATAC_peaks)
#' head(gene_expression)
#' head(TF_footprints)
#' results <- cisREAD(peaks = ATAC_peaks, TFprofile = TF_footprints, genes = c("PRDM1", "MYC", "IRF4", "BCL6")
#' #Look at all candidate coCREs detected by find_coCREs
#' head(results$coCREs)
#' Look at selected (co)CREs
#' head(results$selected_coCREs)
#' Look at significant (co)CREs
#' head(results$significant_coCREs)
#' @export
cisREAD <- function(peaks,
                    TFprofile,
                    RNA,
                    genes,
                    genome = "hg38",
                    distance = 100000,
                    minTFevents=2,
                    coCREcutoff=0.3,
                    weights = c(1,1),
                    lambda = c("lambda_min", "lambda_se"),
                    pval = 0.05,
                    padj = 0.05) {

  print("Mapping candidate cis-regulatory elements to genes")
  peaks2genes <- map_peaks2genes(peaks, genome = genome, distance = distance)

  print("Finding candidate cis-regulatory communities")
  coCRE_list <- lapply(genes,
                       FUN=find_coCREs,
                       peaks2genes = peaks2genes,
                       TFprofile = TFprofile,
                       distance = distance,
                       minTFevents = minTFevents,
                       coCREcutoff = coCREcutoff,
                       weights = c(1,1))

  print("Selecting candidate cis-regulatory communities")
  selected_list <- lapply(coCRE_list,
                          FUN= select_coCREs,
                          RNA = RNA,
                          lambda = lambda,
                          pval = pval,
                          plot = F)

  print("Getting results")
  get_results <- function(runs, list, data) {
    result <- list[[runs]][[data]]
    result
  }

  runs <- 1:length(selected_list)

  all_cCREs <- lapply(runs, FUN = get_results, list = coCRE_list, data = "cCREs")
  all_cCREs <- do.call(rbind, all_cCREs)

  all_coCREs <- lapply(runs, FUN = get_results, list = coCRE_list, data = "coCREs")
  all_coCREs <- do.call(rbind, all_coCREs)

  all_selected <- lapply(runs, FUN = get_results, list = selected_list, data = "selected_coCREs")
  all_selected <- do.call(rbind, all_selected)

  #BH correction
  all_selected$PValue <- replace(all_selected$PValue, is.na(all_selected$PValue), 1)
  gene_pVals <- all_selected %>% dplyr::group_by(Gene) %>% dplyr::summarise(Gene_PValue = min(PValue))
  gene_pVals$Gene_PAdj <- signif(p.adjust(gene_pVals$Gene_PValue, method = "BH"), 3)
  all_selected <- dplyr::left_join(all_selected, gene_pVals, by="Gene")
  all_selected <- all_selected[, -13]
  all_significant <- all_selected[(all_selected$PValue <= pval) & (all_selected$Gene_PAdj <= padj), ]
  rownames(all_significant) <- NULL

  results <- list(cCREs = all_cCREs,
                  coCREs = all_coCREs,
                  selected_coCREs = all_selected,
                  significant_coCREs = all_significant)

  class(results) <- "cisREADResult"
  results
}
