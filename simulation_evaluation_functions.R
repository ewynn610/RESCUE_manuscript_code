## Calculate cell metrics
calc_cell_metrics=function(data, type, meta_dat=NULL, meta_dat_name=NULL){
  ## Calculate library size
  lib_size=colSums(data)
  
  ## Calculate cellwise percentage of 0's
  perc_0=colSums(data==0)/nrow(data)
  
  ## Make dataframe
  df=data.frame(lib_size=lib_size, cell_perc_0=perc_0, data_set=type)
  if(!is.null(meta_dat)){
    df[,meta_dat_name]=meta_dat
  }
  df
}

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

run_icc=function(seurat_obj, time_var, subject_var, 
                 sample_var, baseline_vis,
                 n_genes_choose=800){
  library(dplyr)
  
  ## Rename Meta data
  names(seurat_obj@meta.data)[ names(seurat_obj@meta.data)==time_var]<-"time"
  names(seurat_obj@meta.data)[ names(seurat_obj@meta.data)==subject_var]<-"subject_id"
  names(seurat_obj@meta.data)[ names(seurat_obj@meta.data)==sample_var]<-"sample_number"
  
  ## Extract Counts, meta data
  counts=seurat_obj@assays$RNA$counts
  meta=data.frame(time=factor(seurat_obj$time),subject_id=seurat_obj$subject_id, 
                  sample_number=seurat_obj$sample_number)
  seurat_obj=CreateSeuratObject(counts, meta.data=meta)

  
  ## Filter out genes not detected in a lot of cells
    idx_keep=which(apply(counts, 1, function(x) sum(x>0)>=2))
    seurat_obj=seurat_obj[idx_keep,]
    counts=counts[idx_keep,]

  
  ## Transform data
    seurat_obj=SCTransform(seurat_obj, return.only.var.genes = F, do.correct.umi = F,
                           min_cell=1, method="glmGamPoi")
    vst_sub=seurat_obj@assays$SCT@scale.data
  
  
  rm(seurat_obj)
  
  ## Filter to highly expressed genes
    cpm_counts=edgeR::cpm(counts)
    cpm_counts=cpm_counts[rownames(vst_sub),]
    mean_cpm<-rowMeans(cpm_counts)
    idx_keep=order(mean_cpm, decreasing=T)[1:min(n_genes_choose, length(mean_cpm))]
    vst_sub=vst_sub[idx_keep,]
  
  rm(cpm_counts)
  rm(idx_keep)
  
  ## Run LMM
  print("Running LMM")
  
      res=lmerSeq::lmerSeq.fit(~(1|subject_id)+(1|sample_number), 
                               expr_mat = vst_sub,
                               sample_data=meta)
  
  ## Function to calculate ICC
  calc_icc=function(lmerseq_obj){
    vc=lme4::VarCorr(lmerseq_obj$fit)
    sample_num_var=attr(vc$sample_number, "stddev")^2
    residual_var=sigma(lmerseq_obj$fit)^2
    
      subj_id_var=attr(vc$subject_id, "stddev")^2
      sum_var=sample_num_var+residual_var+subj_id_var
      icc_subj=subj_id_var/sum_var
      icc_samp=sample_num_var/sum_var
    
      data.frame(gene=lmerseq_obj$gene,icc_subj=icc_subj, 
                 icc_samp=icc_samp,
                 re_var_subj=subj_id_var,
                 re_var_samp=sample_num_var,
                 residual_var=residual_var)
    
  }
  
  ## Make dataframe of ICC
  df_icc=lapply(res, function(x){
    calc_icc(x)
  })%>%bind_rows()
  
  
  return(df_icc)
}