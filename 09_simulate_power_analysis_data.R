#!/usr/bin/env -S Rscript --vanilla
#
# Set the log file and job name
#SBATCH -J 09_sim_power_data-%a
#SBATCH -e code/errs/09_sim_power_data-%a
#SBATCH -o code/logs/09_sim_power_data-%a
#SBATCH --array=1-30
#
# Set RAM, nodes, and cores per node
#SBATCH --mem 64G
#SBATCH -n 1
#SBATCH -c 1
#
#######################################
##
## Script name: Power analysis simulation
##
## Purpose of script: Simulate data for power analysis using RAM cells 
##
## Author: Elizabeth Wynn
##
##
#######################################

## load the packages we will need:
library(rescueSim)
library(SingleCellExperiment)

i=as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

#######################################
## 01 - Set simulation parameters
#######################################

## All simulation scenarios
## "biggest" is largest of all three parameters
scenarios=rbind(baseline=c(5, 2, 200),
                subjects_2x=c(10, 2, 200),
                timepoints_2x=c(5, 4, 200),
                cells_2x=c(5, 2, 400),
                biggest=c(10, 4, 400))
colnames(scenarios)=c("nsubjects", "ntimepoints", "ncells")

myScenario=scenarios["biggest",]

## Increase/decrease samp/subj variance
samp_subj_var_scenarios=expand.grid(sim=1:10, samp_subj_var=c("small", "med","large"))
samp_subj_var=samp_subj_var_scenarios$samp_subj_var[i]
samp_subj_var_fac<-ifelse(samp_subj_var=="small", -log(1.5),
                          ifelse(samp_subj_var=="med", 0, log(1.5)))
sim_num=samp_subj_var_scenarios$sim[i]

#######################################
## 02 - Read in simulation parameters
#######################################

## Simulation parameters already computed in 04_simulate_rescueSim_data.R

RAM_params=readRDS("params/mould_RAM_rescueSim_params.RDS")

## Simulate logFC parameters
## I'll draw log2FC's manually so I can use the same values for large/med/small var simulation
## Trying to make them more comparable
set.seed(sim_num^2)
time3<-sample(c(-.35, 0, .35), size = length(getRescueSimParam(RAM_params, "exprsMean")), 
              prob = c(.1, .8, .1), replace=T)
time2<-time3*(2/3)
time1<-time3*(1/3)

## Update parameters with biggest settings and selected samp/subj var
RAM_params <- updateRescueSimParams(
  paramObj = RAM_params,
  paramValues = list(
    nTimepoints = myScenario["ntimepoints"],
    nSubjsPerGroup = myScenario["nsubjects"],
    maxCellsPerSamp = myScenario["ncells"]+100,
    minCellsPerSamp = myScenario["ncells"]-100,
    sampleFacVarMean=RAM_params@sampleFacVarMean+samp_subj_var_fac,
    subjectFacVarMean=RAM_params@subjectFacVarMean+samp_subj_var_fac,
    deLog2FC=list(time1=time1, time2=time2, time3=time3)
  )
)


#######################################
## 03 - Simulate Data
#######################################

## Simulate with biggest settings
mySim <- simRescueData(RAM_params)

## Downsample subjects/timepoints/cells as needed
for(x in 1:4){
  scenario_name=paste0(rownames(scenarios)[x], "_", samp_subj_var, "_samp_subj_var")
  
  samps_per_timepoint_fil=scenarios[x, "nsubjects"]
  timepoints_fil=scenarios[x, "ntimepoints"]
  cells_per_samp_fil=scenarios[x, "ncells"]
  
  print(paste("Samps:", samps_per_timepoint_fil, ", timepoints:",
              timepoints_fil, ", cells: ", cells_per_samp_fil))
  ## Filter to right number of subjs
  subjs_keep=paste0("subject", 1:samps_per_timepoint_fil)
  mySim_new=mySim[,mySim$subjectID %in% subjs_keep]
  
  ## Filter to right number of times
  if(timepoints_fil==2){
    mySim_new$time=as.numeric(gsub("time","", mySim_new$time))
    mySim_new=mySim_new[,mySim_new$time==0| mySim_new$time==3]
    mySim_new$time[mySim_new$time==3]=1
    rowData(mySim_new)=rowData(mySim_new)[,"deLog2FC.time3"]
    names(rowData(mySim_new))="deLog2FC.time1"
  }
  
  ## Filter to right number of cells per sample
  if(cells_per_samp_fil!=scenarios["biggest", "ncells"]){
    max_cells=cells_per_samp_fil+100
    min_cells=cells_per_samp_fil-100
    nsamps=samps_per_timepoint_fil*timepoints_fil
    ncells_per_samp=sample(min_cells:max_cells, nsamps, replace = T)
    names(ncells_per_samp)=unique(mySim_new$sampleID)
    idx_cells_keep=sapply(unique(mySim_new$sampleID), function(x){
      idx_samps=which(mySim_new$sampleID==x)
      ncells=ncells_per_samp[x]
      sample(idx_samps, ncells,replace=F)
    })
    idx_cells_keep=unlist(idx_cells_keep)
    mySim_new=mySim_new[,idx_cells_keep]
  }
  saveRDS(mySim_new, paste0("power_analysis_simulations/RAM_",scenario_name, "_", sim_num, ".RDS" ))
}


