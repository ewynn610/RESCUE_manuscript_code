#!/usr/bin/env -S Rscript --vanilla
#
# Set the log file and job name
#SBATCH -J 05_simulate_splatPop-%a
#SBATCH -e code/errs/err_05_simulate_splatPop-%a
#SBATCH -o code/logs/out_05_simulate_splatPop-%a
#SBATCH --array=1-19
#
# Set RAM, nodes, and cores per node
#SBATCH --mem 512G
#SBATCH -n 1
#SBATCH -c 1
#
#######################################
##
## Script name: Simulate splatPop data
##
## Purpose of script: Simulate data from all papers/celltypes using splatPop
##
## Author: Elizabeth Wynn
##
#######################################


## load up the packages we will need:
library(splatter)
library(rescueSim)
library(dplyr)

## Load custom function to run splatPop
source("code/00_splatPop_custom_function.R")


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

## Read in params to extract batch parameter from
my_rescueSimParams <- readRDS(paste0("params/", my_cell_type, "_rescueSim_params.RDS"))

#######################################
## 02 - Calculate batch parameter
#######################################

## Extract sample batch parameters from rescueSim params
sampleFacVarMean=getRescueSimParam(my_rescueSimParams, "sampleFacVarMean")
sampleFacVarSD<-getRescueSimParam(my_rescueSimParams, "sampleFacVarSD")

## We know v~logNorm(sampleFacVarMean, sampleFacVarSD)
## Then sqrt(v)~lognorm(sampleFacVarMean/2, samplefacVarSD/2)
## E(sqrt(v))=exp(mu + sigma^2/2)=exp(sampleFacVarMean/2 + samplefacVarSD^2/8)
batch.facScale<-exp(sampleFacVarMean/2+sampleFacVarSD^2/8)

#######################################
## 03 - Simulate data with subject variability
#######################################

## Set seed for reproducibility
set.seed(24)
splatpop_dat_subj<-simSplatpop_custom(dat, "sampleID", "subjectID",
                                       non_de_genes = non_de_genes, 
                                       batch.facScale = batch.facScale, 
                                       simSubjCorr = T)

saveRDS(splatpop_dat_subj, paste0("simulations/", my_cell_type, 
                                  "_splatPop_with_subj.RDS"))

rm(splatpop_dat_subj)
gc()

#######################################
## 04 - Simulate data without subject variability
#######################################
## Set seed for reproducibility
set.seed(24)

splatpop_dat_no_subj<-simSplatpop_custom(dat, "sampleID", "subjectID",
                                          non_de_genes = non_de_genes,
                                          simSubjCorr = F)


saveRDS(splatpop_dat_no_subj, paste0("simulations/", my_cell_type,
                                     "_splatPop_no_subj.RDS"))




