#!/usr/bin/env -S Rscript --vanilla
#
#SBATCH -J 07_calc_cell_metrics-%a
#SBATCH -e code/errs/err_07_calc_cell_metrics-%a
#SBATCH -o code/logs/out_07_calc_cell_metrics-%a
#SBATCH --array=1-19
#SBATCH --mem=256G
#SBATCH -c 3

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
library(parallel)
library(SingleCellExperiment)

source("code/00_simulation_metric_functions.R")

#######################################
## 01 - Set-up
#######################################
i <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
file_names <- list.files("data_processed/sce_objs/")
my_file <- file_names[i]
my_cell_type <- sub("_sce.*$", "", my_file)

#######################################
## 02 - Calculate cell metrics
#######################################

# Load empirical data
empirical_sce <- readRDS(paste0("data_processed/sce_objs/", my_cell_type, "_sce.RDS"))
empirical_metrics <- calc_cell_metrics(counts(empirical_sce), "empirical",
                                       data.frame(as.character(empirical_sce$sampleID),
                                                  as.character(empirical_sce$subjectID),
                                                  empirical_sce$time),
                                       c("sampleID","subjectID", "time"))

# Define sims
sim_datasets <- c("rescueSim", "splatPop_no_subj", "splatPop_with_subj")

# Load all sims and compute metrics
all_metrics <- list(empirical_metrics)

rm(empirical_sce, empirical_metrics)
gc()

for (sim_name in sim_datasets) {
  sim_sce <- readRDS(paste0("simulations/", my_cell_type, "_", sim_name, ".RDS"))
  sim_metrics <- calc_cell_metrics(counts(sim_sce), sim_name,
                                   data.frame(as.character(sim_sce$sampleID),
                                              as.character(sim_sce$subjectID),
                                              sim_sce$time),
                                   c("sampleID","subjectID", "time"))
  all_metrics[[length(all_metrics) + 1]] <- sim_metrics
  rm(sim_sce, sim_metrics)
  gc()
}

# Combine all
cell_metrics <- bind_rows(all_metrics)

# Save combined metrics
write.csv(cell_metrics, paste0("simulation_metrics/cell_metrics/", my_cell_type, "_cell_metrics.csv"), row.names = FALSE)

#######################################
## 03 - Univariate KS
#######################################

# Univariate KS
ks_uni <- lapply(c("log_lib_size", "cell_perc_0"), function(x) {
  lapply(sim_datasets, function(y) {
    test <- ks.test(cell_metrics[cell_metrics$data_set == y, x],
                    cell_metrics[cell_metrics$data_set == "empirical", x])
    data.frame(stat = x, data_set = y, ks_stat = test$statistic)
  }) %>% bind_rows()
}) %>% bind_rows()

write.csv(ks_uni, paste0("simulation_metrics/cell_metrics/", my_cell_type, "_cell_metric_ks.csv"), row.names = FALSE)

#######################################
## 04 - Bivariate KS
#######################################

# Bivariate KS: flatten sim + stat combinations
param_grid <- expand.grid(
  stat = list(
    c("log_lib_size", "cell_perc_0")
  ),
  data_set = sim_datasets,
  stringsAsFactors = FALSE
)

ks_bi <- mclapply(seq_len(nrow(param_grid)), function(j) {
  x <- param_grid$stat[[j]]
  y <- param_grid$data_set[j]
  test <- peacock2(
    cell_metrics[cell_metrics$data_set == y, x],
    cell_metrics[cell_metrics$data_set == "empirical", x]
  )
  data.frame(stat = paste(x, collapse = "/"), data_set = y, ks_stat = test)
}, mc.cores = 3) %>% bind_rows()

write.csv(ks_bi, paste0("simulation_metrics/cell_metrics/", my_cell_type, "_cell_metric_ks2.csv"), row.names = FALSE)


