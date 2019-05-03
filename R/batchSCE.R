#' Run batch correction methods on a SCtkExperiment object
#'
#' @param inSCE Input SCtkExperiment object. Required
#' @param batch The name of a column in colData to use as the batch variable
#' @param group The name of a column in colData to use as primary condition variable
#' @param covariates List of other column names in colData to be added to the
#' ComBat model as covariates
#' Combat parameters
#' @param useAssay The assay to use for ComBat. The default is "logcounts"
#' @param par.prior TRUE indicates parametric adjustments will be used, FALSE
#' indicates non-parametric adjustments will be used. Accepted parameters:
#' "Parametric" or "Non-parametric"
#' @param mean.only If TRUE ComBat only corrects the mean of the batch effect
#' @param ref.batch If given, will use the selected batch as a reference for
#' batch adjustment.
#' Combatseq parameters
#' @param shrink If TRUE uses empirical Bayes to estimate batch effect parameters
#' @param gene.subset.n Number of genes to use for empirical Bayes estimation. If not 
#' specified and shrink=TRUE, use 1000 random genes
#' 
#' @return Batch adjusted matrix based on inputs. You can save this matrix into the
#' SCtkExperiment with assay()
#' @export
#' @examples
#' if(requireNamespace("bladderbatch", quietly = TRUE)) {
#'   library(bladderbatch)
#'   data(bladderdata)
#'
#'   #subset for testing
#'   dat <- bladderEset[1:50,]
#'   dat <- as(as(dat, "SummarizedExperiment"), "SCtkExperiment")
#'   mod <- stats::model.matrix(~as.factor(cancer), data = colData(dat))
#'
#'   # parametric adjustment
#'   combat_edata1 <- batchSCE(inSCE = dat, useAssay = "exprs", bea_method='ComBat',
#'                              batch = "batch")
#'   assay(dat, "parametric_combat") <- combat_edata1
#' }
#'
batchSCE <- function(inSCE, batch, group=NULL, covariates=NULL, bea_method,
                     useAssay="logcounts", par.prior="Parametric", mean.only=FALSE, ref.batch=NULL, 
                     shrink=FALSE, gene.subset.n=NULL){
  
  if(bea_method=="ComBat"){
    resassay <- ComBatSCE(inSCE=inSCE, batch=batch, useAssay=useAssay,
                          par.prior=par.prior, covariates=c(group, covariates), 
                          mean.only=mean.only, ref.batch=ref.batch)
  }else if(bea_method=="ComBat-Seq"){
    resassay <- ComBatSeqSCE(inSCE=inSCE, batch=batch, group=group, covariates=covariates, 
                             shrink=shrink, gene.subset.n=gene.subset.n)
  } 
  else{
    stop("Unsupported batch effect adjustment method, ", bea_method)
  }
  
  return(resassay)
}


#' @describeIn batchSCE Perform batch correction with combat
#'
#' @export
#' @examples
#' if(requireNamespace("bladderbatch", quietly = TRUE)) {
#'   library(bladderbatch)
#'   data(bladderdata)
#'
#'   #subset for testing
#'   dat <- bladderEset[1:50,]
#'   dat <- as(as(dat, "SummarizedExperiment"), "SCtkExperiment")
#'   mod <- stats::model.matrix(~as.factor(cancer), data = colData(dat))
#'
#'   # parametric adjustment
#'   combat_edata1 <- ComBatSCE(inSCE = dat, useAssay = "exprs",
#'                              batch = "batch", covariates = NULL)
#'   assay(dat, "parametric_combat") <- combat_edata1
#'
#'   # non-parametric adjustment, mean-only version
#'   combat_edata2 <- ComBatSCE(inSCE = dat, useAssay = "exprs",
#'                              batch = "batch", par.prior = "Non-parametric",
#'                              mean.only = TRUE, covariates = NULL)
#'   assay(dat, "nonparametric_combat_meanonly") <- combat_edata2
#'
#'   # reference-batch version, with covariates
#'   combat_edata3 <- ComBatSCE(inSCE = dat, useAssay = "exprs",
#'                              batch = "batch", covariates = "cancer",
#'                              ref.batch = 3)
#'   assay(dat, "refbatch_combat_wcov") <- combat_edata3
#'   assays(dat)
#' }
#'
ComBatSCE <- function(inSCE, batch, useAssay="logcounts",
                      par.prior="Parametric", covariates=NULL, mean.only=FALSE,
                      ref.batch=NULL){

  #prepare model matrix
  mod <- NULL
  if (length(covariates) > 0){
    mod <- stats::model.matrix(
      stats::as.formula(paste0("~", paste0(covariates, collapse = "+"))),
      data = data.frame(SingleCellExperiment::colData(inSCE)[, covariates,
                                                               drop = FALSE]))
  }

  #prepare parametric
  if (par.prior == "Parametric"){
    par.prior <- TRUE
  } else if (par.prior == "Non-parametric") {
    par.prior <- FALSE
  } else {
    stop("Invalid option given to par.prior. Accepted values are Parametric",
         " and Non-parametric.")
  }

  resassay <-
    sva::ComBat(dat = SummarizedExperiment::assay(inSCE, useAssay),
                batch = SingleCellExperiment::colData(inSCE)[, batch],
                mod = mod, par.prior = par.prior,
                mean.only = mean.only, ref.batch = ref.batch)
  return(resassay)
}


#' @describeIn batchSCE Perform batch correction with combat-seq
#'
#' @export
#' @examples
#' if(requireNamespace("BatchQC", quietly = TRUE)) {
#'   library(BatchQC)
#'   nbatch <- 3
#'   ncond <- 2
#'   npercond <- 5
#'   dat <- rnaseq_sim(ngenes=50, nbatch=nbatch, ncond=ncond, npercond=npercond)
#'   df <- data.frame(batch = rep(1:nbatch, each=ncond*npercond), 
#'                    condition = rep(rep(1:ncond, each=npercond), nbatch))
#'   se <- SummarizedExperiment(assays=list(counts=dat), colData=as(df, "DataFrame"))
#'   dat <- as(se, "SCtkExperiment")
#'
#'   combatseq_counts <- ComBatSeqSCE(inSCE = dat, batch = "batch", group="condition")
#'   assay(dat, "combatseq") <- combatseq_counts
#' }
#'
ComBatSeqSCE <- function(inSCE, batch, group=NULL, covariates=NULL, 
                         shrink=FALSE, gene.subset.n=NULL){
  
  ## Prepare for SCE input
  dat <- SummarizedExperiment::assay(inSCE, 'counts')
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
  mu_hat <- glm_f2$fitted.values
  phi_hat <- do.call(cbind, genewise_disp_lst)
  
  
  ########  In each batch, compute posterior estimation through Monte-Carlo integration  ########  
  if(shrink){
    cat("Apply EB - computing posterior estimates for parameters\n")
    #if(Cpp){mcint_fun <- monte_carlo_int_NB_cpp}else{mcint_fun <- monte_carlo_int_NB}
    mcint_fun <- monte_carlo_int_NB
    monte_carlo_res <- lapply(1:n_batch, function(ii){
      if(ii==1){
        mcres <- mcint_fun(dat=dat[, batches_ind[[ii]]], mu=mu_hat[, batches_ind[[ii]]],
                           gamma=gamma_hat[, ii], phi=phi_hat[, ii], gene.subset.n=gene.subset.n)
      }else{
      invisible(capture.output(mcres <- mcint_fun(dat=dat[, batches_ind[[ii]]], mu=mu_hat[, batches_ind[[ii]]], 
                                                  gamma=gamma_hat[, ii], phi=phi_hat[, ii], gene.subset.n=gene.subset.n)))
      }
      return(mcres)
    })
    names(monte_carlo_res) <- paste0('batch', levels(batch))
    
    gamma_star_mat <- lapply(monte_carlo_res, function(res){res$gamma_star})
    gamma_star_mat <- do.call(cbind, gamma_star_mat)
    phi_star_mat <- lapply(monte_carlo_res, function(res){res$phi_star})
    phi_star_mat <- do.call(cbind, phi_star_mat)
    
    # if(!shrink.disp){
    #   cat("Apply EB shrinkage to mean only\n")
    #   phi_star_mat <- phi_hat
    # }
  }else{
    cat("EB shrinkage off - using GLM estimates for parameters\n")
    gamma_star_mat <- gamma_hat
    phi_star_mat <- phi_hat
  }
  
  
  ########  Obtain adjusted batch-free distribution  ########
  mu_star <- matrix(NA, nrow=nrow(dat), ncol=ncol(dat))
  for(jj in 1:n_batch){
    mu_star[, batches_ind[[jj]]] <- exp(log(mu_hat[, batches_ind[[jj]]])-vec2mat(gamma_star_mat[, jj], n_batches[jj]))
  }
  phi_star <- rowMeans(phi_star_mat)
  
  
  ########  Adjust the data  ########  
  cat("Adjusting the data\n")
  adjust_counts <- matrix(NA, nrow=nrow(dat), ncol=ncol(dat))
  for(kk in 1:n_batch){
    counts_sub <- dat[, batches_ind[[kk]]]
    old_mu <- mu_hat[, batches_ind[[kk]]]
    old_phi <- phi_hat[, kk]
    new_mu <- mu_star[, batches_ind[[kk]]]
    new_phi <- phi_star
    adjust_counts[, batches_ind[[kk]]] <- match_quantiles(counts_sub=counts_sub, 
                                                          old_mu=old_mu, old_phi=old_phi, 
                                                          new_mu=new_mu, new_phi=new_phi)
  }
  
  storage.mode(adjust_counts) <- "integer"
  dimnames(adjust_counts) <- dimnames(dat)
  return(adjust_counts)
}


vec2mat <- function(vec, n_times){
  # expand vector into matrix
  return(matrix(rep(vec, n_times), ncol=n_times, byrow=FALSE))
}


monte_carlo_int_NB <- function(dat, mu, gamma, phi, gene.subset.n){
  # non-parametric monte carlo integration
  pos_res <- lapply(1:nrow(dat), function(i){
    ph <- phi[-i]		
    m <- mu[-i,!is.na(dat[i,])]
    x <- dat[i,!is.na(dat[i,])]
    gamma_sub <- gamma[-i]
    phi_sub <- phi[-i]
    
    # take a subset of genes to do integration - save time
    if(!is.null(gene.subset.n) & is.numeric(gene.subset.n) & length(gene.subset.n)==1){
      if(i==1){cat(sprintf("Using %s random genes for Monte Carlo integration\n", gene.subset.n))}
      mcint_ind <- sample(1:(nrow(dat)-1), gene.subset.n, replace=FALSE)
      m <- m[mcint_ind, ]; ph <- ph[mcint_ind]; gamma_sub <- gamma_sub[mcint_ind]; phi_sub <- phi_sub[mcint_ind]
      G_sub <- gene.subset.n
    }else{
      if(i==1){cat("Using all genes for Monte Carlo integration; the function runs very slow for large number of genes\n")}
      G_sub <- nrow(dat)-1
    }
    
    LH <- sapply(1:G_sub, function(j){prod(dnbinom(x, mu=m[j,], size=1/ph[j]))})
    LH[is.nan(LH)]=0; 
    if(sum(LH)==0 | is.na(sum(LH))){
      return(c(gamma.star=as.numeric(gamma[i]), phi.star=as.numeric(phi[i])))
    }else{
      return(c(gamma.star=sum(gamma_sub*LH)/sum(LH), phi.star=sum(phi_sub*LH)/sum(LH)))
    }
  })
  pos_res <- do.call(rbind, pos_res)
  res <- list(gamma_star=pos_res[, "gamma.star"], phi_star=pos_res[, "phi.star"])	
  return(res)
} 


match_quantiles <- function(counts_sub, old_mu, old_phi, new_mu, new_phi){
  # adjust data by quantile matching
  new_counts_sub <- matrix(NA, nrow=nrow(counts_sub), ncol=ncol(counts_sub))
  for(a in 1:nrow(counts_sub)){
    #print(a)
    for(b in 1:ncol(counts_sub)){
      #print(b)
      if(counts_sub[a, b] <= 1){
        new_counts_sub[a,b] <- counts_sub[a, b]
      }else{
        tmp_p <- pnbinom(counts_sub[a, b]-1, mu=old_mu[a, b], size=1/old_phi[a])
        if(abs(tmp_p-1)<1e-4){
          new_counts_sub[a,b] <- counts_sub[a, b]
          # for outlier count, if p==1, will return Inf values -> use original count instead
        }else{
          new_counts_sub[a,b] <- 1+qnbinom(tmp_p, mu=new_mu[a, b], size=1/new_phi[a])
        }
      }
    }
  }
  return(new_counts_sub)
}

