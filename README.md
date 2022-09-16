# cisREAD
cisREAD (cis-Regualtory Elements Across Differentiation) is an R package designed to identify lineage-specific cis-regulatory elements from ATAC-seq and RNA-seq datasets

# Installation

cisREAD, and its dependancies, can be installed using [devtools](https://devtools.r-lib.org/):

    devtools::install_github("AmberEmmett/cisREAD")
    library(cisREAD)

# Inputs

cisREAD requires the following three input files:

- An RNA-seq count matrix, where rows are genes and columns are samples spanning a cellular lineage, containing log2 normalised counts for Ensembl Gene IDs.

- An ATAC-seq peak file, where the first 6 columns are in BED format and the following columns give log2 normalised ATAC-seq counts in the samples. Samples in RNA-seq and ATAC-seq datasets must match!!! To identify regulatory elements which drive differentiation, we recommend finding differentially-accessible peaks first, using a package like DiffBind or DESeq2.

- A Transcription factor (TF) matrix, where rows are ATAC-seq peaks (the same as in the peak file) and columns are TF occupancy events in a differentiation-stage (e.g. Cell1_TF1, Cell1_TF2, Cell2_TF1, Cell2_TF2 etc.) This is a binary matrix, where 1 indicates the presence of a footprint, for a given TF in a given cell, in the peak, and 0 an absence. In the paper, we derived TF occupancy by discovering enriched transcription factor motifs (with HOMER) and matching these to ATAC-seq footprints (called using HINT-ATAC). Alternate methods for generating the TF matrix are ChIP-seq and motif scans.

# Usage

cisREAD can either be ran on genes one-by-one in a stepwise fashion, or iterated over a large list of genes. We recommend first running the workflow manually on a handful of test genes, to optimise parameters and check that the method works for your data. However if you're in a hurry scroll down to the bottom of the page and take a peak at the last section: 'Applying cisREAD to a list of genes'.

## Step 1. Mapping peaks to genes (map_peaks2genes)

To identify genes near candidate CREs, we first map each ATAC-seq peak to every gene within a user-specified distance, we recommend finding CREs within 100kb of a protein coding TSS. Peak-gene mapping relies on data from ENSEMBL, so please make sure your gene IDs are converted to ENSEMBL identifiers if necessary. Currently only the hg38 annotation of the human genome is supported.

    peaks2genes <- map_peaks2genes(peaks = ATAC_peaks, genome = "hg38", distance = 100000)
    head(peaks2genes)

## Step 2. Finding Communities of Cis-Regulatory Elements (coCREs)

We can now look at the candidate CREs mapped to a gene, and identify those which regulate transcription together. To do this we calculate the chromatin accessibility correlation, and transcription factor occupancy similarity for ever possible pair of candidate CREs, and produce an integrated similarity score. We can then draw edges between candidate CREs whose score exceeds a threshold and construct a network. Infomap community detection, implemented in [igraph](https://igraph.org/r/), is then used to identify communities of connected regulatory elements. These represent co-acting CREs, which are accessible in the same differentiation stages, and occupied by common TFs. To detect coCREs, we use the find_coCREs function.

find_coCREs has a number of adjustable parameters which influence how communities are detected, these are:

- coCREcutoff - The minimum integrated similarity score needed to draw an edge (default 0.3) - lower values result in looser, less similar coCREs, and higher values result in tighter, more similar coCREs.
- minTFevents - The minimum number of TF occupancy events (i.e. 1's in the TF occupancy matrix) needed for each candidate CRE to be taken forward to community detection. This removes candidate CREs from the network which are infrequently TF-bound.
- coCRE_groupings - whether to group CREs based on chromatin accessibility correlation ('ATAC'), transcription factor binding events ('TF'), or integrated similarity scores ('integrated'). The default "integrated" scores are calculated through mutliplying chromatin accessibility correlation by transcription factor binding similarity. 

Lets have a go

    coCREs <- find_coCREs(peaks2genes = peaks2genes, TFprofile = TF_footprints, gene = "PRDM1")
    
  Afterwards we can look at the CREs grouped into communites, the transcription factors bound by each community member, and the accessibility of each community member
  
    #Show all gene-specific candidates
    PRDM1_coCREs$cCREs

    #Show communites of cis-regulatory elements (coCREs)
    PRDM1_coCREs$coCREs
    
    #Show their TF occupancy events
    PRDM1_coCREs$coCRE_TFs

    #Show their chromatin accessibility profiles
    PRDM1_coCREs$coCRE_Accessibility


You can also extract the igraph [communities](https://igraph.org/r/doc/communities.html) and [(plottable) network graph](https://igraph.org/r/doc/plot.common.html)

    PRDM1_coCREs$coCRE_communities
    
    PRDM1_coCRE_network <- PRDM1_coCREs$coCRE_network
    plot(PRDM1_coCRE_network)
    
## Step 3. Selecting gene-specific (co)CREs with LASSO Regression

After identifying CRE communities for each gene, we now narrow down this list by selecting those whose chromatin accessibility best predicts gene expression with LASSO regression. In these models, each candidate (co)CREs chromatin accessibility is an independent variable (X), and each gene's expression is the dependent variable (Y.) LASSO penalises regression co-efficients by assigning them a beta equal to the absolute magnitude of the coefficient. This enables less predictive coefficients to be shrunk towards zero and eliminated from the model, the remaining X variables are considered to be selected.

To determine the degree of shrinkage in a LASSO model we fine-tune the tuning-paramater lambda, to find the model which results in the minimum mean squared cross-validated error 'lambda_min'). If we want a sparser model, we can also use 'lambda se', which is in one standard error of 'lambda_min' With plot == TRUE we can plot the relationship between lambda and number of coefficients. LASSO regression is implemented using [glmnet](https://glmnet.stanford.edu/articles/glmnet.html). 

After model selection by cross validation, and variable selection, we then test the significance of each predictor, accounting for the selection event using a selective inference approach (implemented in [selectiveInference](https://github.com/selective-inference/R-software).

This can all be performed with the select_coCREs function.

    PRDM1_selected <- select_coCREs(coCREs = PRDM1_coCREs, RNA = gene_expression, lambda = 'lambda_min', pval = 0.05)

We can now examine properties of the selected (co)CREs like before.

    #Show selected coCREs
    PRDM1_selected$selected_coCREs

    #Show their TF occupancy events
    PRDM1_selected$selected_coCRE_TFs

    #Show their chromatin accessibility profiles
    PRDM1_selected$selected_coCRE_Accessibility


## Applying cisREAD to a list of genes

To identify cis-regulatory elements across differentiation for a large list of genes, such as genes differentially expressed across the lineage, we can use the all-in-one function 'cisREAD'. The cisREAD function is a wrapper of the above three functions and accepts all paramaters for map_peaks2genes, find_coCREs and select_coCREs (aside from plot, which is switched off). There is also an additional paramater "padj" which determines the FDR cut-off after Benjamini-Hochberg correction for each model created. Use cisREAD to apply the method to a character vector of gene symbols like below:

    #Define a gene list
    gene_list <- c("PRDM1", "XBP1", "IRF4", "PAX5", "BCL6", "SPI1")
    
    #Apply cisREAD
    cisREAD_CREs <- cisREAD(peaks = ATAC_peaks, genome = "hg38", distance = 100000, TFprofile = TF_footprints, genes = gene_list, RNA = gene_expression, lambda = "lambda_min, pval = 0.05, padj = 0.05)
    
  Tada! You can then extract candidate CREs for each gene, candidate CRE communities and selected coCREs.
  
      #Show candidate CREs
      cisREAD_CREs$cCREs
      
      #Show coCREs
      cisREAD_CREs$coCREs
      #Show selected coCREs
      cisREAD_CREs$selected_coCREs
      
      #Show significnat coCREs
      cisREAD_CREs$significant_coCREs
      
