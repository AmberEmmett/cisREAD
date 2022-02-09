#' Find communities of Cis-Regulatory-Elements (coCREs)

#' This function finds communities of co-accessible and co-bound cis-regulatory elements using infomap 
#' community detection

#' @param peaks2genes A dataframe mapping ATAC-seq peaks to nearby genes. Can be created with cisREAD's
#' map_peaks2genes function
#' @param TFprofile A binary matrix where each row is an ATAC-seq peak, and each column indicates the presence (1)
#' or absence (0) of a TF footprint in a given differentiation stage
#' @param gene A string giving a protein coding gene symbol, e.g. "PRDM1" or "MYC"
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
#' @return A find_coCREsResult object, listing candidate CRE communities, their TFs, and accessibilites
#' @example 
#' head(peaks2genes)
#' head(TF_footprints)
#' PRDM1_coCREs <- find_coCREs(peaks2genes = peaks2genes, TFprofile = TF_footprints, gene = "PRDM1", 
#' minTFevents = 2, coCREcutoff = 0.3, weights = c(1,1))
#' #Look at coCREs detected
#' head(PRDM1_coCREs$coCREs)
#' Look at coCRE accessibilities
#' head(PRDM1_coCREs$coCRE_Accessibility)
#' Look at coCRE TF footprints
#' head(PRDM1_coCREs$coCRE_TFs)
#' @export
find_coCREs <- function(peaks2genes, TFprofile, gene, distance=100000, minTFevents=2, coCREcutoff=0.3, weights = c(1,1)) {

  params <- list(gene = gene, distance = distance, minTFevents = minTFevents, coCREcutoff = coCREcutoff)

  gene2peaks <- peaks2genes[peaks2genes$Gene == gene, ]

  if (length(unique(gene2peaks$Ensembl_ID)) > 1) {
    warning(paste0("The gene symbol '", gene, "' matches multiple Ensembl IDs, try a different gene! \n"))
    return(NULL)
  }

  if (nrow(gene2peaks) == 0) {
    warning(paste("uh-oh,", gene, "does not appear in 'peaks2genes'!\n"))
    return(NULL)
  }

  cCREs <- gene2peaks[abs(gene2peaks$Distance) < distance, ]

  if (nrow(cCREs) < 2) {
    warning(paste("Too few candidate CREs for", gene, ", cannot perform community detection. \n"))
    return(NULL)
  }

  row.names(cCREs) <- cCREs$Peak

  #ATAC correlation
  nCells <- length(cCREs) - 9
  cCRESignal <- t(cCREs[, 10:(nCells + 9)])
  cCRECorrelation <- cor(cCRESignal, method="pearson")

  #TF occupancy
  TF <- as.matrix(TFprofile[rownames(TFprofile) %in% colnames(cCRECorrelation), ])
  occupiedcCREs <- names(which(rowSums(TF) >= minTFevents))

  if (length(occupiedcCREs) < minTFevents) {
    warning(paste("Too few occupied candidate CREs for", gene, ", cannot perform community detection.\n"))
    return(NULL)
  }

  TFSimilarity <- as.matrix(1 - arules::dissimilarity(TF, method = "dice"))

  #Network construction
  coCREEdges <- ((weights[1] * TFSimilarity) * (weights[2] * cCRECorrelation))
  diag(coCREEdges) <- 1
  coCREEdges <- coCREEdges[occupiedcCREs, occupiedcCREs]
  coCREEdges[coCREEdges <= coCREcutoff] <- 0
  coCRE_network <- igraph::graph.adjacency(coCREEdges, mode="undirected", weighted=TRUE)

  #Community detection
  coCRE_communities <- igraph::infomap.community(coCRE_network)
  coCRE_groupings <- data.frame(coCRE = coCRE_communities$membership,
                                Peak = coCRE_communities$names)
  coCRE_groupings <- coCRE_groupings[order(coCRE_groupings$coCRE), ]
  rownames(coCRE_groupings) <- NULL

  #Create result dfs
  TFprofile$Peak <- rownames(TFprofile)
  coCRE_TFs <- dplyr::left_join(coCRE_groupings, TFprofile, by = "Peak")
  coCRE_TFs <- coCRE_TFs[, colSums(coCRE_TFs != 0) > 0]

  cCRE_Accessibility <- cCREs[, c(4, 10:(nCells + 9))]
  coCRE_Accessibility <- dplyr::left_join(coCRE_groupings, cCRE_Accessibility, by="Peak")

  rownames(cCREs) <- NULL
  cCREs <- cCREs[, 1:9]

  coCRE_groupings <- dplyr::left_join(coCRE_groupings, cCREs, by = "Peak")
  coCRE_groupings <- coCRE_groupings[, c(3:5, 2, 6:10, 1)]
  results <- list(coCREs = coCRE_groupings,
                  coCRE_TFs = coCRE_TFs,
                  coCRE_Accessibility = coCRE_Accessibility,
                  coCRE_network = coCRE_network,
                  coCRE_communities = coCRE_communities,
                  cCREs = cCREs,
                  params = params)

  class(results) <- "find_coCREsResult"
  results
}
