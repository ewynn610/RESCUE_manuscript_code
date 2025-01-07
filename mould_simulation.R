#######################################
##
## Script name: Simulate Mould Data
##
## Purpose of script: Simulate and evaluate Mould 2020 et al data
##
## Author: Elizabeth Wynn
#######################################

## load up the packages we will need:
library(RESCUE)
library(Seurat)
library(Matrix.utils)
library(SingleCellExperiment)
library(edgeR)
library(dplyr)
library(patchwork)


myPath="/home/wynne/dissertation/"
source(paste0(myPath,"simulation/code/00_simulation_evaluation_functions.R"))

cell_types <- c("RecAM", "T cell", "Cycling", "RAM")

for(i in 1:length(cell_types)){
  #######################################
  ## 01 - Read in data
  #######################################
  my_cell_type=cell_types[i]
  print(my_cell_type)
  
  ## Read in baseline data
  dat_baseline=readRDS(paste0("mould_data/mould_", my_cell_type, "_baseline.RDS"))
  dat_baseline=as.SingleCellExperiment(dat_baseline)
  
  ## Read in paired data (only non-DEGs included)
  paired_dat=readRDS(paste0("mould_data/mould_", my_cell_type, "paired_non_degs.RDS"))
  
  #######################################
  ## 02 - Estimate params
  #######################################

  ##############Estimate params besides multiplicative factor parameters using baseline data #################

  ## Set params equal to 1 so they won't get estimated
  myParams=RescueParams(sampleFacVarSD=1, sampleFacVarMean=1,subjectFacVarMean=1, subjectFacVarSD=1,maxCellsPerSamp=1,
                        minCellsPerSamp=1)
  myParams <- estRescueParams(
    sce = dat_baseline,
    paramObj=myParams,
    sampleVariable = "sample_id",
    subjectVariable = "subject_id",
    timepointVariable = "time_bin",
    groupVariable = NULL,
    cellParamsByCondition = T
  )
  
  ############## Estimate multiplicative factor params using paired data #################
  
   myParams@sampleFacVarSD=myParams@sampleFacVarMean=myParams@subjectFacVarMean=
    myParams@subjectFacVarSD=myParams@maxCellsPerSamp=myParams@minCellsPerSamp=numeric(0)
  myParams_full <- estRescueParams(
    sce = paired_dat,
    paramObj=myParams,
    sampleVariable = "sample_id",
    subjectVariable = "subject_id",
    timepointVariable = "time_bin",
    groupVariable = NULL, 
    cellParamsByCondition = T
  )
  
   saveRDS(myParams_full,paste0(myPath, "mould_simulations/mould_", my_cell_type, "_params.RDS"))
  
  #######################################
  ## 03 - Simulate data
  #######################################

  mySim <- simRescueData(myParams_full)

  saveRDS(mySim,paste0(myPath,"mould_simulations/mould_",
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
  
  write.csv(cell_metrics,paste0("mould_simulation_summaries/mould_",
                                my_cell_type, "_cell_metrics.csv"))
  
  ############## Gene level metrics #################
  gene_metrics=calc_gene_metrics(counts(mySim), "simulation")
  gene_met$gene=rownames(gene_met)
  gene_met
  
  
  write.csv(gene_metrics,paste0("mould_simulation_summaries/mould_",
                                my_cell_type, "_gene_metrics.csv"))
  
  ############## Calculate icc #################
  icc_dat=run_icc(mySim, 
                  time_var = "time_bin", subject_var = "subjectID",
                  sample_var = "sampleID", baseline_vis="time0")
  
  write.csv(icc_dat,paste0("mould_simulation_summaries/mould_", 
                           my_cell_type, "_icc.csv"))
  
}
