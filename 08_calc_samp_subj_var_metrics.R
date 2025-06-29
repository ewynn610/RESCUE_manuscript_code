#!/usr/bin/env -S Rscript --vanilla
#
#SBATCH -J 08_calc_samp_subj_var_metrics-%a
#SBATCH -e code/errs/err_08_calc_samp_subj_var_metrics-%a
#SBATCH -o code/logs/out_08_calc_samp_subj_var_metrics-%a
#SBATCH --array=1-19
#SBATCH --mem=512G
#SBATCH -c 16

#######################################
## Script name: Cell-level benchmarking statistics
##
## Purpose of script: Compute cell-level stats + KS tests
##
## Author: Elizabeth Wynn
#######################################

library(Seurat)
library(Peacock.test)
library(dplyr)
library(ggplot2)
library(SingleCellExperiment)

source("code/00_simulation_metric_functions.R")

## Set for scTrans
options(future.globals.maxSize = 8000 * 1024^2)

#######################################
## 01 - Read in data
#######################################

i <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
file_names <- list.files("data_processed/sce_objs/")
my_file <- file_names[i]
my_cell_type <- sub("_sce.*$", "", my_file)

## Mapping empirical timepoints to time0, time1, time2
## Note - mapping alphabetically so they match ncells in RESCUE - not temporally
empirical_time0=c(khoo="Acute", lambo="Diagnosis", mould="baseline")
empirical_time1 <- c(knoo="Convalescent", lambo="Relapse", "follow-up")
empirical_time2 <- c(lambo="Remission")

## Function for loading data
load_data_ls <- function(my_cell_type) {
  
  ## Read in empirical data
  empirical = readRDS(paste0("data_processed/sce_objs/", my_cell_type, "_sce.RDS"))
  
  ## Switch time variable so it looks like simulations
  empirical$time=case_when(
    empirical$time %in% empirical_time0 ~ "time0",
    empirical$time %in% empirical_time1 ~ "time1",
    empirical$time %in% empirical_time2 ~ "time2",
    TRUE ~ NA_character_
  )
  
  ## Read in Rescue data data
  rescue_sim = readRDS(paste0("simulations/", my_cell_type, "_rescueSim.RDS"))
  
  ## Subset empirical and rescue to empirical non-de genes
  non_de_genes = readRDS(paste0("data_processed/non_de_genes/", my_cell_type, "_non_de_genes.RDS"))
  empirical <- empirical[non_de_genes,]
  rescue_sim <- rescue_sim[non_de_genes,]
  
  ## Read in splatPop data
  ## Match splatPop genes to non_de genes with empirical means
  ## Filter to matched genes
  splatPop_dat<-lapply(c("splatPop_with_subj", "splatPop_no_subj"), function(sim){
    dat = readRDS(paste0("simulations/", my_cell_type, "_", sim, ".RDS"))
    
    ## Take lognorm of empirical non-de genes, splatPop sim genes
    emp_dat_rowMeans=rowMeans(logcounts(scater::logNormCounts(empirical[non_de_genes,])))
    my_dat_rowMeans=rowMeans(logcounts(scater::logNormCounts(dat)))
    
    non_de_match=c()
    for (emp_mean in emp_dat_rowMeans) {
      match_gene=names(which.min(abs(my_dat_rowMeans - emp_mean)))
      non_de_match=c(non_de_match,match_gene)
      my_dat_rowMeans=my_dat_rowMeans[names(my_dat_rowMeans)!=match_gene]
    }
    
    ## Need at least a count in each cell
    while(any(colSums(counts(dat)[non_de_match,])==0)){
      cells0<-which(colSums(counts(dat)[non_de_match,])==0)
      genes_keep=unique(sapply(cells0, function(y){
        rownames(dat)[which(counts(dat[,y])!=0)][1]
      }))
      non_de_match=c(genes_keep,non_de_match[-1:-length(genes_keep)])
    }
    
    dat <- dat[non_de_match,]
  })
  names(splatPop_dat)<-c("splatPop_with_subj", "splatPop_no_subj")
  list(
    empirical=empirical,
    splatPop_with_subj = splatPop_dat$splatPop_with_subj,
    splatPop_no_subj = splatPop_dat$splatPop_no_subj,
    rescue_sim = rescue_sim
    
  )
}

## Load data
data_ls <- load_data_ls(my_cell_type)

#######################################
## num - Set up color palette for plot
#######################################


## Pallete for Tsne
my_pal <- c(
  "#8dd3c7", "#4d6f5e", "#1a3c3a",  
  "#ffffb3", "#e6cc00", "#b3a600",  
  "#bebada", "#8b6e9c", "#532d72",  
  "#fb8072", "#d84c38", "#9d2320", 
  "#80b1d3", "#3e7a98", "#1e446e",
  "#fdb462", "#e68a2f", "#b76d24",
  "#b3de69", "#4c9c2f", "#246b15",
  "#fccde5", "#f174af", "#c80078",
  "#d9d9d9", "#9e9e9e", "#606060",
  "#bc80bd", "#8a408c", "#5a1f67",
  "#ccebc5", "#4d7f5d", "#1e4e34",
  "#ffed6f", "#e6b500", "#cc8e00",
  "#e41a1c", "#b31512", "#820e0a",
  "#377eb8", "#236093", "#0f4c70", 
  "#4daf4a", "#368f2e", "#1e6420", 
  "#ff7f00", "#e56e00", "#b35b00",
  "#984ea3", "#702785", "#4c1f5f"  
)

## Naming according to biggest number of subjects/timepoints
names(my_pal)=with(expand.grid(time = 1:3, subj = 1:17),
                   paste0("Subj. ", subj, ", Time ", time))

#######################################
## num - Calculate all metrics for each dataset
#######################################

tsne_ls<-list()
sill_width_ls <- list()
icc_ls <- list()
cms_ls <- list()
ks_ls <- list()
for(name_dat in names(data_ls)){
  dat<-data_ls[[name_dat]]
  
  ## Remove genes with all 0's
  dat <- dat[rowMeans(counts(dat))!=0,]
  
  ## Create Seurat obj
  dat_seurat=CreateSeuratObject(counts(dat), meta.data = data.frame(colData(dat)))
  dat_seurat<-SCTransform(dat_seurat, min_cell=0)
  dat_seurat<-RunPCA(dat_seurat)
  
  ## Run Tsne
  tsne_ls[[name_dat]]<-make_tsne(dat_seurat, type=name_dat, col_pal = my_pal)
  
  ## Run sil width
  sill_width_ls[[name_dat]]<-run_sill_width(dat_seurat, name_dat)
  
  rm(dat_seurat)
  gc()
  
  ## Run scTransform but need to include ALL cells
  normcounts(dat)=Seurat::SCTransform(counts(dat),cell.attr = colData(dat), 
                                      return.only.var.genes = F, do.correct.umi = F,
                                      min_cell=0)$y
  
  icc_ls[[name_dat]] <- run_icc(dat, name_dat)
  
  cms_ls[[name_dat]] <- run_cms(dat, name_dat)
  
  if(name_dat!="empirical"){
    cms_ks<-lapply(c("sample_cms", "subject_cms"), function(x){
      test=ks.test(cms_ls[[name_dat]][,x], 
                   cms_ls[["empirical"]][,x])
      data.frame(stat=x, data_set=name_dat, ks_stat=test$statistic)
    })%>%bind_rows()
    
    sw_ks<-lapply(c("sw_subj", "sw_samp"), function(x){
      test=ks.test(sill_width_ls[[name_dat]][,x], 
                   sill_width_ls[["empirical"]][,x])
      data.frame(stat=x, data_set=name_dat, ks_stat=test$statistic)
    })%>%bind_rows()
    
    icc_ks<-lapply(c("icc_subj", "icc_samp"), function(x){
      test=ks.test(icc_ls[[name_dat]][,x], 
                   icc_ls[["empirical"]][,x])
      data.frame(stat=x, data_set=name_dat, ks_stat=test$statistic)
    })%>%bind_rows()
    
    ks_ls[[name_dat]]=bind_rows(list(cms_ks, sw_ks, icc_ks))
  }
  
  
}

sill_width_df <- bind_rows(sill_width_ls)
icc_df <- bind_rows(icc_ls)
cms_df <- bind_rows(cms_ls)
ks_df <- bind_rows(ks_ls)

write.csv(sill_width_df, paste0("simulation_metrics/samp_subj_var_metrics/", my_cell_type, "_sw.csv"), row.names = FALSE)
write.csv(icc_df, paste0("simulation_metrics/samp_subj_var_metrics/", my_cell_type, "_icc.csv"), row.names = FALSE)
write.csv(cms_df, paste0("simulation_metrics/samp_subj_var_metrics/", my_cell_type, "_cms.csv"), row.names = FALSE)
write.csv(ks_df, paste0("simulation_metrics/samp_subj_var_metrics/", my_cell_type, "_subj_samp_var_metric_ks.csv"), row.names = FALSE)


saveRDS(tsne_ls, paste0("simulation_metrics/samp_subj_var_metrics/", my_cell_type, "_tsne_list.RDS"))









