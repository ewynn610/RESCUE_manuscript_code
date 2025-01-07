#######################################
##
## Script name: Khoo simulation
##
## Purpose of script: Simulate data using Khoo 2023 dataset as reference
##
## Author: Elizabeth Wynn
######################################
## load up the packages we will need:
library(RESCUE)
library(Matrix.utils)
library(Seurat)
library(SingleCellExperiment)
library(edgeR)
library(dplyr)


#######################################
## Simulate each cell type of interest
#######################################

cell_types=c("B", "CD4 T", "CD8 T", "NK")

for(i in 1:length(cell_types)){
  #######################################
  ## 01 - Read in data
  #######################################
  my_cell_type=cell_types[i]
  print(my_cell_type)
  
  ## Data was formatted in format_khoo_data.R
  dat=readRDS(paste0("khoo_data/khoo_", my_cell_type, 
                     ".RDS"))
  
  dat=as.SingleCellExperiment(dat)
  
  
  #######################################
  ## 02 - Find non-de genes
  #######################################
  set.seed(24)
  counts=counts(dat)
  
  counts_agg=t(aggregate.Matrix(t(counts), dat$sample, fun="sum"))
  
  meta_data=data.frame(sample_id=dat$sample, subject_id=dat$Patient, 
                       time_bin=dat$Status)
  meta_data_agg=meta_data[!duplicated(meta_data$sample_id),]
  counts_agg=counts_agg[,order(colnames(counts_agg))]
  meta_data_agg=meta_data_agg[order(meta_data_agg$sample_id),]
  identical(colnames(counts_agg), meta_data_agg$sample_id)
  
  ############## Run edgeR with subject as fixed effect #################
  
  design_mat=model.matrix(~time_bin+subject_id, data=meta_data_agg)
  dge <- edgeR::DGEList(counts = counts_agg)
  dge<-edgeR::calcNormFactors(dge)
  dge <- edgeR::estimateDisp(dge, design_mat)
  
  ############## Hypothesis testing for time effect #################
  fit_edgeR <- edgeR::glmFit(dge, design_mat)
  my_coef=paste0("time_bin", levels(factor(meta_data_agg$time_bin))[2])
  res_tab<-edgeR::glmLRT(fit_edgeR, coef=my_coef)$table
  res_tab$padj=p.adjust(res_tab$PValue, method = "BH")
  
  ############## Save if gene is non_de or not #################
  genes_keep_edgeR=rownames(res_tab)[res_tab$padj>.05 &abs(res_tab$logFC)<.1]
  
  
  #######################################
  ## 03 - Estimate params
  #######################################
  
  myParams <- estRescueParams(
    sce = dat,
    sampleVariable = "sample",
    subjectVariable = "Patient",
    timepointVariable = "Status",
    groupVariable = NULL, nonDEGs = genes_keep_edgeR,
    cellParamsByCondition = T
  )
  
  saveRDS(myParams,paste0("khoo_simulations/khoo_", my_cell_type, "_params.RDS"))
  
  #######################################
  ## 04 - Simulate data
  #######################################
  
  mySim <- simRescueData(myParams)
  
  saveRDS(mySim,paste0("khoo_simulations/khoo_",
                       my_cell_type, "_sim.RDS"))
  
  #######################################
  ## 05 - Compute summary statistics
  #######################################
  
  ############## cell level metrics #################
  cell_metrics=calc_cell_metrics(counts(mySim), "simulation",
                                 data.frame(as.character(mySim$sampleID),
                                            as.character(mySim$subjectID),
                                            mySim$time),
                                 c("sampleID","subjectID", "time"))
  
  write.csv(cell_metrics,paste0("khoo_simulation_summaries/khoo_",
                                my_cell_type, "_cell_metrics.csv"))
  
  ############## Gene level metrics #################
  gene_metrics=calc_gene_metrics(counts(mySim), "simulation")
  gene_met$gene=rownames(gene_met)
  gene_met
  
  
  write.csv(gene_metrics,paste0("khoo_simulation_summaries/khoo_",
                                my_cell_type, "_gene_metrics.csv"))
  
  ############## Calculate icc #################
  icc_dat=run_icc(mySim, 
                  time_var = "time_bin", subject_var = "subjectID",
                  sample_var = "sampleID", baseline_vis="time0")
  
  write.csv(icc_dat,paste0("khoo_simulation_summaries/khoo_", 
                           my_cell_type, "_icc.csv"))
  
  
}
