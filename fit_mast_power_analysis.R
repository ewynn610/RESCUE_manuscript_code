#######################################
## Script name: MAST models power simulation
##
## Purpose of script: Fit MAST models for power simulation 
##
## Author: Elizabeth Wynn
#######################################

## load the packages we will need:
library(MAST)

#######################################
## num - Set up scenario conditions
#######################################
scenarios=rbind(base=c(5, 2, 200),
                subjects=c(10, 2, 200),
                timepoints=c(5, 4, 200),
                cells=c(5, 2, 400))
colnames(scenarios)=c("nsubjects", "ntimepoints", "ncells")

for(i in length(scenarios)){
  scenario_name=rownames(scenarios)[i]
  scenario_name=ifelse(scenario_name=="baseline", scenario_name, paste0("double_", scenario_name))
  myScenario=scenarios[i,]
  
  nsims=10
  
  
  
  #######################################
  ## num - Loop through data
  #######################################
  
  ############## Set options #################
  perc0_fil=.9
  scale_fac=1000000
  cores=16
  
  for (j in 1:10) {
    print(j)
    sim=readRDS(paste0("power_analysis_simulations/RAM_",scenario_name, "_", i, ".RDS" ))
    
    sim=sim[apply(counts(sim), 1, function(x) sum(x==0)/length(x))<=perc0_fil,]
    sim$time=as.numeric(gsub("time", "", sim$time))
    
    de=rowData(sim)
    normcounts(sim)<-apply(counts(sim), 2, function(x){
      log2((scale_fac*x/sum(x))+1)
    })
    
    sca<- FromMatrix(normcounts(sim), colData(sim))
    cdr <-colSums(SummarizedExperiment::assay(sca)>0)
    colData(sca)$cngeneson <- scale(cdr)
    
    options(mc.cores=cores)
    zlmCond <- MAST::zlm(form=~ cngeneson + time + (1 | sampleID)+(1|subjectID),
                         sca, method='glmer',ebayes = F,
                         strictConvergence = T,
                         parallel = T)
    
    
    options(mc.cores=cores)
    raw_res=suppressMessages(MAST::summary(zlmCond, doLRT="time", parallel = T))
    my_sum_tab_raw=data.frame(raw_res$datatable)
    my_sum_tab_raw=my_sum_tab_raw[my_sum_tab_raw$component=="H",]
    my_sum_tab=data.frame(p_val_raw=my_sum_tab_raw$Pr..Chisq.)
    rownames(my_sum_tab)=my_sum_tab_raw$primerid
    
    my_sum_tab$p_val_adj=p.adjust(my_sum_tab$p_val_raw, method = "BH")
    my_sum_tab=data.frame(merge(my_sum_tab, de, by=0))
    my_sum_tab$de=ifelse(my_sum_tab$deLogFC==0, F, T)
    
    saveRDS(zlmCond, paste0("power_analysis_results/RAM_",scenario_name, "_", "_fits.RDS" ))
    saveRDS(my_sum_tab, paste0("power_analysis_results/RAM_",scenario_name, "_", j, "_sum.RDS" ))
    
  }
}

