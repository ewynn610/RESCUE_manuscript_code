#!/usr/bin/env -S Rscript --vanilla
#
# Set the log file and job name
#SBATCH -J 10_run_mast_power_analysis-%a
#SBATCH -e code/errs/10_run_mast_power_analysis-%a
#SBATCH -o code/logs/10_run_mast_power_analysis-%a
#SBATCH --array=1-12
#
# Set RAM, nodes, and cores per node
#SBATCH --mem 96G
#SBATCH -n 1
#SBATCH -c 16
#
#######################################
##
## Script name: MAST power analysis
##
## Purpose of script: Run MAST on simulated data for power analysis
##
## Author: Elizabeth Wynn
##
#######################################


## load the packages we will need:
library(MAST)
library(SingleCellExperiment)

#######################################
## 01 - Set up scenario conditions
#######################################
i=as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

## Scenarios
scenarios=c("baseline", "subjects_2x", "timepoints_2x", "cells_2x")
samp_subj_var=c("small", "med", "large")
scenario_mat<-expand.grid(scenarios=scenarios, samp_subj_var=samp_subj_var)

## File name
file_name=paste(scenario_mat$scenarios[i], scenario_mat$samp_subj_var[i], "samp_subj_var", sep="_")

#######################################
## 02 - Run Mast on all simulations
#######################################

############## Set options #################
perc0_fil=.9
scale_fac=1000000
cores=16
nsims=10

## Loop through all simulations
for (j in 1:nsims) {
  sim=readRDS(paste0("power_analysis_simulations/RAM_",file_name, "_", j, ".RDS" ))

  ## Remove genes with no counts in >90% of cells  
  sim=sim[apply(counts(sim), 1, function(x) sum(x==0)/length(x))<=perc0_fil,]
  sim$time=as.numeric(gsub("time", "", sim$time))
  
  ## Save whether DE or not
  de=rowData(sim)
  de$de=ifelse(de[,1]==0, F, T)

  ## Perform log2(cpm+1) normalization
  normcounts(sim)<-apply(counts(sim), 2, function(x){
    log2((scale_fac*x/sum(x))+1)
  })
  
  ## Set up data for mast
  sca<- FromMatrix(as.matrix(normcounts(sim)), colData(sim))
  
  ## Calculate cellular detection rate
  cdr <-colSums(SummarizedExperiment::assay(sca)>0)
  colData(sca)$cngeneson <- scale(cdr)
  
  ## Fit models with cdr and time fixed effects, sampleID, subjectID random effects
  options(mc.cores=cores)
  zlmCond <- MAST::zlm(form=~ cngeneson + time + (1 | sampleID)+(1|subjectID),
                       sca, method='glmer',ebayes = F,
                       strictConvergence = T,
                       parallel = T)

  ## Test for time effect
  options(mc.cores=cores)
  raw_res=suppressMessages(MAST::summary(zlmCond, doLRT="time", parallel = T))
  
  ## Format results - filter to hurdle results
  my_sum_tab_raw=data.frame(raw_res$datatable)
  my_sum_tab_raw=my_sum_tab_raw[my_sum_tab_raw$component=="H",]
  my_sum_tab=data.frame(p_val_raw=my_sum_tab_raw$Pr..Chisq.)
  rownames(my_sum_tab)=my_sum_tab_raw$primerid
  
  ## Calculate adjusted p-value
  my_sum_tab$p_val_adj=p.adjust(my_sum_tab$p_val_raw, method = "BH")
  
  ## Merge with DE
  my_sum_tab=data.frame(merge(my_sum_tab, de, by=0))
  
  ## Save
  saveRDS(my_sum_tab, paste0("power_analysis_summaries/RAM_",file_name, "_", j, "_sum.RDS" ))
}
