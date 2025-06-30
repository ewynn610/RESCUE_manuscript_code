#######################################
##
## Script name: Caluclate power metrics
##
## Purpose of script: Calculate power/FDR/T1E from simulations in power analysis
##
## Author: Elizabeth Wynn
##
#######################################
## load the packages we will need:
library(dplyr)
library(stringr)

#######################################
## 01 - Calculate power/fdr/t1e for all simulations
#######################################

all_files <- list.files("power_analysis_summaries", 
                        full.names = TRUE, pattern = "sum")



power_summary=lapply(all_files, function(x){
  
  ## Read in data
  summ=readRDS(x)
  
  ## Sim number
  sim_num <- str_extract(x, "(?<=_)(\\d+)(?=_sum.RDS)")
  
  ## Calculate power, fdr, t1e at 0.05 threhshold
  pthresh = 0.05
  power=sum(summ$p_val_adj<pthresh&summ$de, na.rm = T)/sum(summ$de)
  fdr=sum(summ$p_val_adj<pthresh&!summ$de, na.rm = T)/sum(summ$p_val_adj<pthresh, na.rm=T)
  t1e=sum(summ$p_val_raw<pthresh&!summ$de, na.rm = T)/sum(!summ$de)
  data.frame(power=power, fdr=fdr,t1e=t1e, sig.level=pthresh, sim_num=sim_num)%>%
    ## Extract scenario and between samp/subj variance setting
    mutate(scenario=case_when(grepl("baseline",x)~"Baseline Case",
                             grepl("cells_2x", x)~"2x # of Cells",
                             grepl("timepoints_2x", x)~"2x # of Timepoints",
                             grepl("subjects_2x", x)~ "2x # of Subjects",
                             grepl("cells_4x", x)~"4x # of Cells"),
           samp_subj_var=case_when(grepl("small_samp_subj_var", x)~"small_var",
                                   grepl("large_samp_subj_var", x)~"large_var",
                                   grepl("med_samp_subj_var", x) ~"med_var" ))
})%>%bind_rows()

write.csv(power_summary, "power_analysis_summaries/power_analysis_metrics.csv",
          row.names = F)
