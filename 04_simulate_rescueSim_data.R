#!/usr/bin/env -S Rscript --vanilla
#
# Set the log file and job name
#SBATCH -J 04_simulate_rescueSim-%a
#SBATCH -e code/errs/err_04_simulate_rescueSim-%a
#SBATCH -o code/logs/out_04_simulate_rescueSim-%a
#SBATCH --array=1-19
#
# Set RAM, nodes, and cores per node
#SBATCH --mem 256G
#SBATCH -c 1
#
#######################################
##
## Script name: Simulate rescueSim data
##
## Purpose of script: Simulate data from all papers/celltypes using rescueSim
##
## Author: Elizabeth Wynn
##
#######################################


## load up the packages we will need:
library(rescueSim)

#######################################
## 01 - Read in data for celltype
#######################################

## Running a job for each celltype
i=as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

## Cell type data saved in:
##    01_format_mould_data.R, 02_format_khoo_data.R, 03_format_lambo_data.R
file_names=list.files("data_processed/sce_objs/")

## Separate job for each file
my_file<- file_names[i]
my_cell_type <- sub("_sce.*$", "", my_file)

## Read in sce and non DE genes
dat <- readRDS(paste0("data_processed/sce_objs/", my_cell_type, "_sce.RDS"))
non_de_genes <- readRDS(paste0("data_processed/non_de_genes/", my_cell_type, 
                               "_non_de_genes.RDS"))

# #######################################
# ## 02 - Estimate params 
# #######################################

## Set seed for reproducibility
set.seed(24)

## Estimate params - use non_de_genes for batch var params
myParams <- estRescueSimParams(
  sce = dat,
  sampleVariable = "sampleID",
  subjectVariable = "subjectID",
  timepointVariable = "time",
  groupVariable = NULL,nonDEGs = non_de_genes,
  cellParamsByCondition = T
)

## Save param values
saveRDS(myParams,paste0("params/", my_cell_type, "_rescueSim_params.RDS"))

gc()

#######################################
## 03 - Simulate data
#######################################

## Simulate data
mySim <- simRescueData(myParams)

## Save results
saveRDS(mySim,paste0("simulations/", my_cell_type, "_rescueSim.RDS"))



