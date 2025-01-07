#######################################
##
## Script name: Simulate data for power analysis
##
## Purpose of script: 
##
## Author: Elizabeth Wynn
##
##
#######################################

## load the packages we will need:
library(RESCUE)

#######################################
## 01 - Set simulation parameters
#######################################

## All simulation scenarios
## "biggest" is largest of all three parameters
scenarios=rbind(base=c(5, 2, 200),
      subjects=c(10, 2, 200),
      timepoints=c(5, 4, 200),
      cells=c(5, 2, 400),
      biggest=c(10, 4, 400))
colnames(scenarios)=c("nsubjects", "ntimepoints", "ncells")

myScenario=scenarios["biggest",]

nsims=10


#######################################
## 02 - Read in simulation parameters
#######################################

## Simulation parameters already computed in home/dissertation/simulation/code/simulation_redo/03_simulate_kara.R

RAM_params=readRDS(paste0(myPath, "simulation/params/kara_RAM_params.RDS"))

## 
RAM_params <- updateRescueParams(
  paramObj = RAM_params,
  paramValues = list(
    nTimepoints = myScenario["ntimepoints"],
    nSubjsPerGroup = myScenario["nsubjects"],
    maxCellsPerSamp = myScenario["ncells"]+100,
    minCellsPerSamp = myScenario["ncells"]-100,
    propDE=.2,
    deLogFC=.35
  )
)


for (i in 1:nsims) {
  print(i)
  set.seed(i)
  mySim <- simRescueData(RAM_params)
  
  for(x in 1:4){
    scenario_name=rownames(scenarios)[x]
    scenario_name=ifelse(scenario_name=="baseline", scenario_name, paste0("double_", scenario_name))
    
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
    saveRDS(mySim_new, paste0("power_analysis_simulations/RAM_",scenario_name, "_", i, ".RDS" ))
  }
  
  
}
