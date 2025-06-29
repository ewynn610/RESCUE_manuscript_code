#######################################
##
## Script name: 
##
## Purpose of script: 
##
## Author: Elizabeth Wynn
##
## Date Created: 2025-06-28
##
#######################################
##
## Notes:
##   
##
#######################################

## load the packages we will need:


#######################################
## num - Define function for calculating gene metrics
#######################################

calc_gene_metrics=function(data, type){
  ## Log normalize data
  norm_dat=Seurat::NormalizeData(data, scale.factor=1000000)
  
    ## Calculate genewise percentage 0's
    perc_0=rowSums(data==0)/ncol(data)
    
    ## Log normalized means
    avg_log_cpm=norm=rowMeans(as.matrix(norm_dat))
    
    ## Log normalized variance
    var_log_cpm=apply(as.matrix(norm_dat), 1, var)
    
    ## Save dataframe of all metrics
    data.frame( gene_perc_0=perc_0, avg_log_cpm=avg_log_cpm,
                var_log_cpm=var_log_cpm,data_set=type)
  
}

#######################################
## 02 - Define function for calculating cell metrics
#######################################

calc_cell_metrics=function(data, type, meta_dat=NULL, meta_dat_name=NULL){
  ## Calculate library size
  log_lib_size=log(colSums(data))
  
  
  ## Calculate cellwise percentage of 0's
  perc_0=colSums(data==0)/nrow(data)
  
  ## Make dataframe
  df=data.frame(log_lib_size=log_lib_size, cell_perc_0=perc_0, data_set=type)
  if(!is.null(meta_dat)){
    df[,meta_dat_name]=meta_dat
  }
  df
}

#######################################
## 03 - Make tsne plots
#######################################
make_tsne <- function(seurat_obj, type, group_var = "sample_id", dims = 30, col_pal = NULL) {
  # Map subject IDs to clean labels
  subj_unique <- unique(seurat_obj$subjectID)
  subj_id_map <- setNames(paste0("Subj. ", seq_along(subj_unique)), subj_unique)
  seurat_obj@meta.data$subject_id <- subj_id_map[as.character(seurat_obj$subjectID)]
  
  # Map time to friendly labels
  time_map <- c("time0" = "Time 1", "time1" = "Time 2", "time2" = "Time 3")
  seurat_obj@meta.data$time <- time_map[seurat_obj$time]
  
  # Combine subject_id and time into sample_id
  seurat_obj$sample_id <- paste0(seurat_obj$subject_id, ", ", seurat_obj$time)
  
  # Run t-SNE
  seurat_obj <- RunTSNE(seurat_obj, dims = 1:dims, check_duplicates = FALSE)
  
  # Build plot
  tsne_plot <- DimPlot(seurat_obj, group.by = group_var, reduction = "tsne") +
    ggtitle(type)
  
  if (!is.null(col_pal)) {
    tsne_plot <- tsne_plot + scale_color_manual(values = col_pal)
  }
  
  return(tsne_plot)
}


#######################################
## 04 - Calculate sillohette width
#######################################

run_sill_width<-function(seurat_obj, type, seed=24){
  ## Sil width doesn't work well if to big...if more than 30,000 cells, downsample
  if(ncol(seurat_obj)>30000){
  set.seed(seed)
  idx<-sample(1:ncol(seurat_obj), size = 30000, replace = F)
  seurat_obj <- seurat_obj[,idx]
  }
  
  ## Dist matrix
  mtx <- dist(seurat_obj[["pca"]]@cell.embeddings)
  
  ## Subject and sample silhouette width
  res_subj <- cluster::silhouette(as.integer(factor(seurat_obj$subjectID)), mtx)
  res_samp <- cluster::silhouette(as.integer(factor(seurat_obj$sampleID)), mtx)
  
  ## Return summary
  data.frame(sw_subj=res_subj[,"sil_width"], sw_samp=res_samp[,"sil_width"],
             data_set=type, sampleID=seurat_obj$sampleID,
             subjectID=as.character(seurat_obj$subjectID), time=seurat_obj$time)
}

#######################################
## 05 - Run ICC
#######################################

run_icc=function(sce,type, n_genes_choose=800){
 
  ## Filter to highly expressed genes
  vst<-normcounts(sce)
  mean_norm<-rowMeans(vst)
  idx_keep=order(mean_norm, decreasing=T)[1:min(n_genes_choose, length(mean_norm))]
  vst_sub=vst[idx_keep,]
  ## Run LMM
  print("Running LMM")
  
  res=lmerSeq::lmerSeq.fit(~(1|subjectID)+(1|sampleID), 
                           expr_mat = vst_sub,
                           sample_data=data.frame(colData(sce)), parallel = T,
                           cores=16)
  
  ## Function to calculate ICC
  calc_icc=function(lmerseq_obj){
    vc=lme4::VarCorr(lmerseq_obj$fit)
    sample_num_var=attr(vc$sampleID, "stddev")^2
    residual_var=sigma(lmerseq_obj$fit)^2
    
    subj_id_var=attr(vc$subjectID, "stddev")^2
    sum_var=sample_num_var+residual_var+subj_id_var
    icc_subj=subj_id_var/sum_var
    icc_samp=sample_num_var/sum_var
    
    data.frame(gene=lmerseq_obj$gene,icc_subj=icc_subj, 
               icc_samp=icc_samp,
               re_var_subj=subj_id_var,
               re_var_samp=sample_num_var,
               residual_var=residual_var, 
               type=type)
    
  }
  
  ## Make dataframe of ICC
  df_icc=lapply(res, function(x){
    calc_icc(x)
  })%>%bind_rows()
  
  
  return(df_icc)
}

#######################################
## 06 - Run CMS
#######################################


run_cms<-function(dat, type){
  
  ## Calculate k based on empirical data
  n_empirical=table(data_ls$empirical$sampleID)
  k=max(10, round(median(n_empirical[n_empirical != 0])))
  
  ## Sample cms
  sample_cms<-CellMixS::cms(dat,k=k, k_min=10, group="sampleID", assay_name="normcounts", 
                  BPPARAM = BiocParallel::MulticoreParam(16))
  
  ## Subject cms
  n_empirical=table(data_ls$empirical$subjectID)
  k=max(10, round(median(n_empirical[n_empirical != 0])))
  subject_cms<-CellMixS::cms(dat,k=k, k_min=10, group="subjectID", assay_name="normcounts", 
                   BPPARAM =BiocParallel::MulticoreParam(16))
  
  data.frame(sample_cms=sample_cms$cms, subject_cms=subject_cms$cms,
             data_set=type, sampleID=dat$sampleID,
             subjectID=as.character(dat$subjectID), time=dat$time)
  
}



