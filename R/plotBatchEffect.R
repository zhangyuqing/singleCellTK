#' Visualize batch effect with different kinds of diagnostic plots
#'
#' @param inSCE Input SCtkExperiment object. Required
#' @param useAssay Indicate which assay to use for PCA. Default is "logcounts"
#' @param batch The column in the annotation data that corresponds to batch.
#' Required
#' @param condition The column in the annotation data that corresponds to
#' condition. Optional
#' @param log_dat Default FALSE, if TRUE log transform the assay
#' @param method Type of plot to generate
#' 
#' @return A plot of degree of batch effect.
#' @export
#' @examples
#' if(requireNamespace("bladderbatch", quietly = TRUE)) {
#'   library(bladderbatch)
#'   data(bladderdata)
#'   dat <- as(as(bladderEset, "SummarizedExperiment"), "SCtkExperiment")
#'   plotBatchEffect(dat, useAssay="exprs", batch="batch", condition = "cancer",
#'                   method="Explained Variation")
#' }
plotBatchEffect <- function(inSCE, useAssay="logcounts", batch,
                            condition=NULL, covariates=NULL, log_dat=FALSE, method){
  if(method=="Explained Variation"){
    a <- plotBatchVariance(inSCE, useAssay=useAssay, batch=batch, 
                           condition=condition, log_dat=log_dat)
    return(a)
  }else if(method %in% c("Mean Batch Effect Estimates", "Dispersion Batch Effect Estimates")){
    batch_params <- estimateParams(inSCE, useAssay=useAssay, batch=batch, 
                                   group=condition, covariates=covariates)
    if(!is.null(batch_params)){
      if(method == "Mean Batch Effect Estimates"){
        batch_params_plt <- batch_params$gamma_hat
      }else if(method == "Dispersion Batch Effect Estimates"){
        batch_params_plt <- batch_params$phi_hat
      }
      colnames(batch_params_plt) <- gsub("batch", "", colnames(batch_params_plt))
      batch_params_plt_mlt <- melt(batch_params_plt, varnames=c("Genes", "Batch"))[,-1]
      a <- ggplot(batch_params_plt_mlt, aes(x=Batch, y=value)) +
        geom_violin() +
        labs(y=method)
      return(a)
    }
  }
  return(NULL)
}


plotBatchVariance <- function(inSCE, useAssay="logcounts", batch,
                              condition=NULL, log_dat=FALSE){
  # Plot the percent of the variation that is explained by batch and condition in the data
  nlb <- nlevels(as.factor(SingleCellExperiment::colData(inSCE)[, batch]))
  if (nlb <= 1){
    batchMod <- matrix(rep(1, ncol(inSCE)), ncol = 1)
  } else {
    batchMod <- stats::model.matrix(
      ~as.factor(SingleCellExperiment::colData(inSCE)[, batch]))
  }
  if (is.null(condition)){
    stop("condition required for now")
  } else {
    nlc <- nlevels(as.factor(
      SingleCellExperiment::colData(inSCE)[, condition]))
    if (nlc <= 1){
      condMod <- matrix(rep(1, ncol(inSCE)), ncol = 1)
    } else {
      condMod <- stats::model.matrix(
        ~as.factor(SingleCellExperiment::colData(inSCE)[, condition]))
    }
  }

  mod <- cbind(condMod, batchMod[, -1])
  
  dat <- SummarizedExperiment::assay(inSCE, useAssay)
  if(log_dat){
    print("logging data")
    dat <- log(dat + 0.25)
  }
  condTest <- batchqc_f.pvalue(dat, mod, batchMod)
  batchTest <- batchqc_f.pvalue(dat, mod, condMod)

  r2Full <- condTest$r2Full
  condR2 <- batchTest$r2Reduced
  batchR2 <- condTest$r2Reduced
  explainedVariation <- round(cbind(`Full (Condition+Batch)` = r2Full,
                                     Condition = condR2,
                                     Batch = batchR2), 5) * 100
  exVarM <- reshape2::melt(explainedVariation)
  colnames(exVarM) <- c("Gene", "Model", "Percent.Explained.Variation")
  exVarM$Model <- factor(exVarM$Model)
  a <- ggplot2::ggplot(exVarM, ggplot2::aes_string("Model", "Percent.Explained.Variation")) +
    ggplot2::geom_violin(ggplot2::aes_string(fill = "Model")) +
    ggplot2::geom_boxplot(width = .1) +
    ggplot2::xlab("Model") +
    ggplot2::ylab("Percent Explained Variation") +
    ggplot2::scale_fill_manual(values = RColorBrewer::brewer.pal(9, "Set1"),
                               guide = FALSE)
  return(a)
}


batchqc_f.pvalue <- function(dat, mod, mod0) {
  # F-test (full/reduced model) and returns R2 values
  # (full/reduced) as well.
  mod00 <- matrix(rep(1, ncol(dat)), ncol = 1)
  n <- dim(dat)[2]
  m <- dim(dat)[1]
  df1 <- dim(mod)[2]
  df0 <- dim(mod0)[2]
  p <- rep(0, m)
  
  resid <- dat - dat %*% mod %*% solve(t(mod) %*% mod) %*% t(mod)
  rss1 <- rowSums(resid * resid)
  rm(resid)
  
  resid0 <- dat - dat %*% mod0 %*% solve(t(mod0) %*% mod0) %*% t(mod0)
  rss0 <- rowSums(resid0 * resid0)
  rm(resid0)
  
  resid00 <- dat - dat %*% mod00 %*% solve(t(mod00) %*% mod00) %*% t(mod00)
  rss00 <- rowSums(resid00 * resid00)
  rm(resid00)
  
  r2Full <- 1 - rss1 / rss00
  r2Reduced <- 1 - rss0 / rss00
  
  p <- 1
  if (df1 > df0)  {
    fstats <- ((rss0 - rss1) / (df1 - df0)) / (rss1 / (n - df1))
    p <- 1 - stats::pf(fstats, df1 = (df1 - df0), df2 = (n - df1))
  }
  return(list(p = p, r2Full = r2Full, r2Reduced = r2Reduced))
}


estimateParams <- function(inSCE, batch, useAssay, group=NULL, covariates=NULL, 
                           shrink=FALSE, gene.subset.n=NULL){
  ## Prepare for SCE input
  dat <- SummarizedExperiment::assay(inSCE, useAssay)
  if(!is.integer(dat[1,1])){return(NULL)}
  batch <- SingleCellExperiment::colData(inSCE)[, batch]
  if(is.null(group)){
    full_mod <- FALSE
  }else{
    full_mod <- TRUE
    group <- SingleCellExperiment::colData(inSCE)[, group]
  }
  covar_mod <- NULL
  if(length(covariates)>0){
    covar_mod <- stats::model.matrix(
      stats::as.formula(paste0("~", paste0(covariates, collapse = "+"))),
      data = data.frame(SingleCellExperiment::colData(inSCE)[, covariates,
                                                             drop = FALSE]))
  }
  if(shrink & is.null(gene.subset.n)){gene.subset.n <- min(nrow(dat)-1, 1000)}
  
  dge_obj <- edgeR::DGEList(counts=dat)
  ## Prepare characteristics on batches
  batch <- as.factor(batch)
  n_batch <- nlevels(batch)  # number of batches
  batches_ind <- lapply(1:n_batch, function(i){which(batch==levels(batch)[i])}) # list of samples in each batch  
  n_batches <- sapply(batches_ind, length)
  #if(any(n_batches==1)){mean_only=TRUE; cat("Note: one batch has only one sample, setting mean.only=TRUE\n")}
  n_sample <- sum(n_batches)
  cat("Found",n_batch,'batches\n')
  
  ## Make design matrix 
  # batch
  batchmod <- model.matrix(~-1+batch)  # colnames: levels(batch)
  # covariate
  group <- as.factor(group)
  if(full_mod & nlevels(group)>1){
    cat("Using full model in ComBat-seq.\n")
    mod <- model.matrix(~group)
  }else{
    cat("Using null model in ComBat-seq.\n")
    mod <- model.matrix(~1, data=as.data.frame(t(dat)))
  }
  # drop intercept in covariate model
  if(!is.null(covar_mod)){covar_mod <- covar_mod[, !apply(covar_mod, 2, function(x){all(x==1)})]}
  # bind with biological condition of interest
  mod <- cbind(mod, covar_mod)
  # combine
  design <- cbind(batchmod, mod)
  
  ## Check for intercept in covariates, and drop if present
  check <- apply(design, 2, function(x) all(x == 1))
  #if(!is.null(ref)){check[ref]=FALSE} ## except don't throw away the reference batch indicator
  design <- as.matrix(design[,!check])
  cat("Adjusting for",ncol(design)-ncol(batchmod),'covariate(s) or covariate level(s)\n')
  
  ## Check if the design is confounded
  if(qr(design)$rank<ncol(design)){
    #if(ncol(design)<=(n_batch)){stop("Batch variables are redundant! Remove one or more of the batch variables so they are no longer confounded")}
    if(ncol(design)==(n_batch+1)){stop("The covariate is confounded with batch! Remove the covariate and rerun ComBat")}
    if(ncol(design)>(n_batch+1)){
      if((qr(design[,-c(1:n_batch)])$rank<ncol(design[,-c(1:n_batch)]))){stop('The covariates are confounded! Please remove one or more of the covariates so the design is not confounded')
      }else{stop("At least one covariate is confounded with batch! Please remove confounded covariates and rerun ComBat")}}
  }
  
  ## Check for missing values in count matrix
  NAs = any(is.na(dat))
  #if(NAs){cat(c('Found',sum(is.na(dat)),'Missing Data Values\n'),sep=' ')}
  
  
  ########  Estimate gene-wise dispersions within each batch  ########
  cat("Estimating dispersions\n")
  ## Estimate common dispersion within each batch as an initial value
  disp_common <- sapply(1:n_batch, function(i){
    if(n_batches[i]==1){
      stop("ComBat-seq doesn't support 1 sample per batch yet!")
    }else if((n_batches[i] <= ncol(design)-ncol(batchmod)+1) | qr(mod[batches_ind[[i]], ])$rank < ncol(mod)){ 
      # not enough residual degree of freedom
      return(edgeR::estimateGLMCommonDisp(dat[, batches_ind[[i]]], design=NULL, subset=nrow(dat)))
      #as.matrix(design[batches_ind[[i]], (n_batch+1):ncol(design)]),
    }else{
      return(edgeR::estimateGLMCommonDisp(dat[, batches_ind[[i]]], design=mod[batches_ind[[i]], ], subset=nrow(dat)))
    }
  })
  
  ## Estimate gene-wise dispersion within each batch 
  genewise_disp_lst <- lapply(1:n_batch, function(j){
    if(n_batches[j]==1){
      stop("ComBat-seq doesn't support 1 sample per batch yet!")
    }else if((n_batches[j] <= ncol(design)-ncol(batchmod)+1) | qr(mod[batches_ind[[j]], ])$rank < ncol(mod)){
      # not enough residual degrees of freedom - use the common dispersion
      return(rep(disp_common[j], nrow(dat)))
    }else{
      return(edgeR::estimateGLMTagwiseDisp(dat[, batches_ind[[j]]], design=mod[batches_ind[[j]], ], 
                                           dispersion=disp_common[j], prior.df=0))
    }
  })
  names(genewise_disp_lst) <- paste0('batch', levels(batch))
  
  ## construct dispersion matrix
  phi_matrix <- matrix(NA, nrow=nrow(dat), ncol=ncol(dat))
  for(k in 1:n_batch){
    phi_matrix[, batches_ind[[k]]] <- vec2mat(genewise_disp_lst[[k]], n_batches[k]) #matrix(rep(genewise_disp_lst[[k]], n_batches[k]), ncol=n_batches[k])
  }
  
  
  ########  Estimate parameters from NB GLM  ########
  cat("Fitting the GLM model\n")
  glm_f <- edgeR::glmFit(dge_obj, design=design, dispersion=phi_matrix, prior.count=1e-4) #no intercept - nonEstimable; compute offset (library sizes) within function
  alpha_g <- glm_f$coefficients[, 1:n_batch] %*% as.matrix(n_batches/n_sample) #compute intercept as batch-size-weighted average from batches
  new_offset <- t(vec2mat(edgeR::getOffset(dge_obj), nrow(dat))) +   # original offset - sample (library) size
    vec2mat(alpha_g, ncol(dat))  # new offset - gene background expression
  # getOffset(dge_obj) is the same as log(dge_obj$samples$lib.size)
  glm_f2 <- edgeR::glmFit.default(dge_obj$counts, design=design, dispersion=phi_matrix, offset=new_offset, prior.count=1e-4) 
  
  #beta_hat <- glm_f2$coefficients[, (n_batch+1):ncol(design)]
  gamma_hat <- glm_f2$coefficients[, 1:n_batch]
  phi_hat <- do.call(cbind, genewise_disp_lst)
  
  return(list(gamma_hat=gamma_hat, phi_hat=phi_hat))
}

vec2mat <- function(vec, n_times){
  # expand vector into matrix
  return(matrix(rep(vec, n_times), ncol=n_times, byrow=FALSE))
}

