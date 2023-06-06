#' Selects cis-regulatory elements and cis-regulatory communities to regulate a gene

#' This function finds CREs and coCREs whose chromatin accessibility best predicts the expresion of a target
#' gene using LASSO regression. This should be used downstream of 'find_coCREs'

#' @param coCREs A find_coCREsResult object created with cisREAD's find_coCREs function
#' @param RNA A gene expression matrix where each row is a gene's ENSEMBL ID (if not then please convert first)
#' and each column gives the gene's expression (log2 normalised counts) in a given sample
#' @param lambda The value of lambda used for model selection during cross-validation, one of 
#' "lambda_min" (default) and "lambda_se". "lambda_min" will construct the model with the mean
#' minimum cross-validated squared error, "lambda_se" will construct the model within 1 standard error 
#' of "lambda_min", and will result in a sparser model
#' @param pval The p-value cutoff for assigning statistical significance. coCRE predictors with
#' p < this value will be considered statistically significant
#' @param plot Indicates whether to construct glmnet plots during LASSO cross-validation, will visualise the 
#' relationship between lambda and the number of co-efficients selected (default = FALSE)
#' @return A select_coCREsResult object, listing selected CRE communities, their TFs and accessibilities, alongside
#' the glmnet model object, inputs and coefficients
#' 
#'
#' @example 
#' head(gene_expression)
#' #First perform find_coCREs
#' PRDM1_coCREs <- find_coCREs(peaks2genes = peaks2genes, TFprofile = TF_footprints, gene = "PRDM1", 
#' minTFevents = 2, coCREcutoff = 0.3, weights = c(1,1))
#' #Next select coCREs which predict gene expression
#' PRDM1_selected_coCREs <- select_coCREs(coCREs = PRDM1_coCREs, RNA = gene_expression, lambda = "lambda_min", pval = 0.05, plot = F)
#' #Look at significant CREs
#' head(PRDM1_selected_coCREs$selected_coCREs)
#' @export
select_coCREs <- function(coCREs, RNA, lambda = c("lambda_min", "lambda_se") , pval = 0.05, plot = FALSE) {

  lambda <- rlang::arg_match(lambda)

  if(is.null(coCREs)) {
    return(NULL)
  }

  set.seed(1)

  params <- list(gene = coCREs$params$gene, distance = coCREs$params$distance, lambda, pval)
  nCells <- ncol(coCREs$coCRE_Accessibility) - 2

  #Standardised x matrix
  coCRESignal <- aggregate(coCREs$coCRE_Accessibility[-(1:2)], list(coCRE = coCREs$coCRE_Accessibility$coCRE), mean)
  rownames(coCRESignal) <- coCRESignal$coCRE
  z_coCRESignal <- scale(t(coCRESignal[, -1]))

  #Standardised y  vector
  gene_exp <- RNA[unique(coCREs$cCREs$Ensembl_ID), ]
  gene_exp <- na.omit(gene_exp)

  if(nrow(gene_exp) == 0) {
    warning(paste0(coCREs$params$gene, " is not in your gene expression matrix!\n"))
    return(NULL)
  }

  z_gene_exp <-  as.vector(t(scale(t(gene_exp))))

  #LASSO model
  model <- tryCatch(

    expr = {
      model <- suppressWarnings(glmnet::glmnet(z_coCRESignal, z_gene_exp, standardize = FALSE))
    },

    error = function(model) {
      warning(paste0("Cannot create a LASSO model for ", coCREs$params$gene, ".\n"))
      warning(" Here's the error message from glmnet: ")
      warning(model)
      return(NA)
    }
  )

  if (is.na(head(model, 1))) return(NULL)

  #LASSO CV
  cv_model <- tryCatch(
    expr = {
      cv_model <- suppressWarnings(glmnet::cv.glmnet(z_coCRESignal,
                                                     z_gene_exp,
                                                     standardize = FALSE,
                                                     nfolds = 10,
                                                     nlambda = 100,
                                                     thresh = 1e-20,
                                                     grouped = FALSE))
    },

    error = function(cv_model) {
      warning(paste0("Cannot cross validate the model for ", coCREs$params$gene, ".\n"))
      warning(" Here's the error message from glmnet: ")
      warning(cv_model)
      return(NA)
    }
  )

  if (is.na(head(cv_model, 1))) return(NULL)

  #Plot
  if (plot == TRUE) {
    par(mfrow = c(2, 1))
    plot(model, label = TRUE, xvar="lambda")
    plot(cv_model, xlab = "Log Lambda")
  }

  #Final model
  opt_lambda <- ifelse(lambda == "lambda_min", cv_model$lambda.min, cv_model$lambda.1se)
  final_model <- suppressWarnings(glmnet::glmnet(z_coCRESignal,
                                                 z_gene_exp,
                                                 standardize = FALSE,
                                                 lambda = opt_lambda,
                                                 thresh = 1e-20))

  coef <- glmnet::coef.glmnet(final_model, exact = TRUE, standardize = FALSE, thresh = 1e-20)
  selected_coCREs <- data.frame(coCRE = as.double(rownames(coef)[-1]),
                                Coefficient = signif(as.vector(coef)[-1], 3))
  selected_coCREs <- selected_coCREs[selected_coCREs$Coefficient != 0, ]

  if (nrow(selected_coCREs) == 0 ) {
    warning(paste0("No candidate CREs have been selected for ", coCREs$params$gene,".\n"))
    return(NULL)
  }

  #Selective inference
  beta <- glmnet::coef.glmnet(final_model,
                              x = z_coCRESignal,
                              y = z_gene_exp,
                              exact = TRUE,
                              standardize = FALSE,
                              thresh = 1e-20)[-1]

  sigma <- suppressWarnings(selectiveInference::estimateSigma(z_coCRESignal,z_gene_exp, standardize=FALSE))

  pval <- tryCatch(
    expr = {
      pval <- suppressWarnings(selectiveInference::fixedLassoInf(x = z_coCRESignal,
                                                                 y = z_gene_exp,
                                                                 beta = beta,
                                                                 lambda = opt_lambda*nCells,
                                                                 alpha = pval,
                                                                 sigma = sigma$sigmahat))
    },

    error = function(pval) {
      warning(paste0("Cannot calculate significance for CREs selected for ", coCREs$params$gene, ". Output selected CREs will have p-values of NA.\n"))
      warning("Here's the error message from selectiveInference: ")
      warning(pval)
      return(NA)
    }
  )

  suppressWarnings(
    if (is.null(pval)) {

      selected_coCREs$PValue <- NA
      selected_coCREs <- dplyr::left_join(selected_coCREs, coCREs$coCREs, by = "coCRE")

    } else if (is.na(head(pval, 1))) {

      selected_coCREs$PValue <- NA
      selected_coCREs <- dplyr::left_join(selected_coCREs, coCREs$coCREs, by = "coCRE")

    } else {

      coCRE_pval <- data.frame(coCRE = as.double(pval$vars), PValue = signif(pval$pv, 3))
      selected_coCREs <- dplyr::left_join(selected_coCREs, coCRE_pval, by = "coCRE")
      selected_coCREs <- dplyr::left_join(selected_coCREs, coCREs$coCREs, by = "coCRE")
    }
  )

  #Create result dfs
  selected_coCREs <- selected_coCREs[c(4:12, 1:3)]

  selected_coCRE_TFs <- coCREs$coCRE_TFs[coCREs$coCRE_TFs$Peak %in% selected_coCREs$Peak, ]
  selected_coCRE_TFs <- selected_coCRE_TFs[, colSums(selected_coCRE_TFs != 0) > 0]
  rownames(selected_coCRE_TFs) <- NULL

  selected_coCRE_Accessibility <- coCREs$coCRE_Accessibility[coCREs$coCRE_Accessibility$Peak %in% selected_coCREs$Peak, ]
  rownames(selected_coCRE_Accessibility) <- NULL


  results <- list(selected_coCREs = selected_coCREs,
                  selected_coCRE_TFs = selected_coCRE_TFs,
                  selected_coCRE_Accessibility = selected_coCRE_Accessibility,
                  x_input = z_coCRESignal,
                  y_input = z_gene_exp,
                  final_model = final_model,
                  model_coefficients = coef,
                  pval = pval,
                  params = params)

  class(results) <- "selectCREsResult"
  results
}
