#######################################
##
## Script name: Custom splatPop function
##
## Purpose of script: Write function to simulate and 
##
## Author: Elizabeth Wynn
##
#######################################
##
## Notes: We will simulate splatpop using just sample level correlation
##        and also using batches to get sample and subject level correlation.
##
#######################################

## load the packages we will need:

#######################################
## 01 - Define main function
#######################################

simSplatpop_custom <- function(dat,
                               sampleID,
                               subjectID,
                               non_de_genes,
                               batch.facScale,
                               simSubjCorr) {
  
  ## Get data information
  counts_full <- counts(dat)
  col_dat <- colData(dat)
  
  nsamps <- length(unique(col_dat[[sampleID]]))
  nsubjs <- length(unique(col_dat[[subjectID]]))
  
  ## If simulating subject correlation, aggregate over subject, else sample
  agg_var <- if (simSubjCorr) subjectID else sampleID
  
  ## Calculate Pseudo means across agg_var
  pseudo_means <- counts(scuttle::aggregateAcrossCells(
    dat[non_de_genes, ], 
    ids = col_dat[[agg_var]], 
    statistics = "mean"
  ))
  
  ## Select biggest sample for parameter estimation
  ## only use one sample to get right variability estimates
    samp_tab <- table(col_dat[[sampleID]])
    simSample <- names(which.max(samp_tab))
  dat_sub <- dat[, col_dat[[sampleID]] == simSample]
  counts_one_samp <- counts(dat_sub)
  
  ## Estimate splatPop parameters
  params_est <- estimate_splatpop_params(counts_one_samp, pseudo_means)
  
  ## Estimate cell number parameters for gamma distribution
  n_cells <- as.numeric(table(col_dat[[sampleID]]))
  ncells_fit <- fitdistrplus::fitdist(n_cells, "gamma", method = "mle", lower = c(0, 0))
  
  ## Set splatPop parameters
    ## no eqtl, 
    ## draw from ncells using gamma params, 
    ## pop.quant.norm=F because we used mean aggregation
  params_est <- setParams(params_est, list(
    eqtl.n = 0,
    nCells.shape = ncells_fit$estimate["shape"],
    nCells.rate = ncells_fit$estimate["rate"],
    pop.quant.norm = FALSE
  ))
  
  ## Run if simulating subject level correlation
  if (simSubjCorr) {
    ## Draw the # of cells per sample (batch) from a gama using gamma estimates
    batch_cells <- pmax(1, round(rgamma(nsamps, ncells_fit$estimate["shape"], ncells_fit$estimate["rate"])))
    
    ## Set params
      ## nCells.sample = F (we include cells per sample in batch_cells)
      ## batch.size = 1, limit batch to one subject ("population")
      ## batch.facLoc = 0: batch factors centered around log(0)=1
      ## batch facScale: use batch.facScale derived from rescueSim
    params_est <- setParams(params_est, list(
      nCells.sample = FALSE,
      batchCells = batch_cells,
      batch.size = 1,
      batch.facLoc = 0,
      batch.facScale = batch.facScale
    ))
    
    ## How many "populations" for vcf
    n_samples_vcf <- nsubjs
  } else {
    
    ## Set params: estimate ncells per sample
    params_est <- setParams(params_est, list(
      nCells.sample = TRUE
    ))
    ## How many "populations" for vcf
    n_samples_vcf <- nsamps
  }
  
  ## Simulate
  splatpop_data <- splatPopSimulate(params_est, vcf = mockVCF(n.samples = n_samples_vcf))
  colnames(splatpop_data) <- paste0("cell_", seq_len(ncol(splatpop_data)))
  
  ## Build metadata
  splatpop_data <- add_metadata(splatpop_data, simSubjCorr, nsubjs, nsamps, sampleID)
  
  ## Get rid of unneccessary pieces (save on size)
  splatpop_data@assays@data$BatchCellMeans<-splatpop_data@assays@data$BaseCellMeans<-
    splatpop_data@assays@data$BCV<-splatpop_data@assays@data$CellMeans<-
    splatpop_data@assays@data$TrueCounts<-NULL
  
  return(splatpop_data)
}

#######################################
## 02 - Helper functions
#######################################

estimate_splatpop_params <- function(counts_one_samp, pseudo_means) {
  
  ## Sometimes get error if pop.cv.bins is too high
  ## Go lower until there's no error
  pop_cv_bins <- 10
  while (pop_cv_bins > 0) {
    params <- newSplatPopParams(pop.cv.bins = pop_cv_bins)
    message("Trying pop.cv.bins = ", pop_cv_bins)
    try_res <- tryCatch({
      splatPopEstimate(counts = as.matrix(counts_one_samp), 
                       means = as.matrix(pseudo_means), 
                       params = params)
    }, error = function(e) {
      message("Error: ", e$message)
      NULL
    })
    if (!is.null(try_res) && !any(is.na(try_res@pop.cv.param))) {
      return(try_res)
    }
    pop_cv_bins <- pop_cv_bins - 1
  }
  stop("splatPopEstimate failed for all pop.cv.bins values.")
}

add_metadata <- function(splatpop_data, simSubjCorr, nsubjs, nsamps, sampleID) {
  if (simSubjCorr) {
    splatpop_data$subjectID <- gsub("sample_0|sample_", "subject", splatpop_data$Sample)
    metdat <- data.frame(
      cell_type = my_cell_type,
      sampleID = gsub("Batch", "sample", splatpop_data$Batch),
      subjectID = splatpop_data$subjectID,
      row.names = colnames(splatpop_data)
    )
    metdat_single <- metdat[!duplicated(metdat$sampleID), ] %>%
      dplyr::group_by(subjectID) %>%
      dplyr::mutate(time = paste0("time", dplyr::row_number() - 1)) %>%
      dplyr::ungroup() %>%
      as.data.frame()
    rownames(metdat_single) <- metdat_single$sampleID
    metdat$time <- metdat_single[metdat$sampleID, "time"]
  } else {
    splatpop_data$sampleID <- gsub("_0|_", "", splatpop_data$Sample)
    metdat <- data.frame(
      cell_type = my_cell_type,
      sampleID = splatpop_data$sampleID,
      cell = colnames(splatpop_data)
    )
    samps <- unique(metdat$sampleID)
    subjs <- rep(paste0("subj_", seq_len(nsubjs)), length.out = length(samps))
    time_bin <- rep(paste0("time", seq_len(ceiling(length(samps)/nsubjs)) - 1), each = nsubjs)[seq_along(samps)]
    subj_metdat <- data.frame(sampleID = samps, subjectID = subjs, time = time_bin)
    metdat <- tibble::column_to_rownames(merge(metdat, subj_metdat), "cell")
  }
  colData(splatpop_data) <- S4Vectors::DataFrame(metdat[colnames(splatpop_data), ])
  return(splatpop_data)
}
