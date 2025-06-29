#######################################
##
## Script name: Save Khoo data
##
## Purpose of script: Format and save Khoo data for simulation
##
## Author: Elizabeth Wynn
##
#######################################

## load up the packages we will need:
library(Seurat)
library(edgeR)
library(scuttle)


#######################################
## 01 - Load in data
#######################################

## Seurat object downloaded from geo database (GSE196456)
load("data_raw/khoo_raw.Rdata")

## Filter to only children
dat=Merged_Seuobj_l1[,Merged_Seuobj_l1$Age_group=="Kid"]

## Save sample identifier
dat$sampleID <- paste(dat$Patient,dat$Status, sep="_")

## Save variables with new names for cosistency
dat$subjectID <- dat$Patient
dat$time <- dat$Status

## Create new seurat object to get rid of junk
dat=CreateSeuratObject(counts=dat@assays$RNA@counts,
                       meta.data=dat[[]])

rm(Merged_Seuobj_l1)
gc()

#######################################
## 02 - Save object for each cell types
#######################################

cell_types=unique(dat$Final_anno.l1)

## Remove "others" category
cell_types=cell_types[cell_types !="others"]

## Loop through cell types
for(my_cell_type in cell_types){
  ## Subset to single cell type
  print(my_cell_type)
  dat_fil=subset(dat, Final_anno.l1==my_cell_type)
  
  ## Want at least 5 samples with 25 cells in each sample
  if(sum(rowSums(table(dat_fil$subjectID,dat_fil$time)>=25)==2)<5){
    next
  }
  
  ## Remove genes that are 0's for all cells
  dat_fil=dat_fil[rowSums(dat_fil@assays$RNA$counts)!=0,]
  
  dat_fil=as.SingleCellExperiment(dat_fil)
  
  saveRDS(dat_fil, paste0("data_processed/sce_objs/khoo_", my_cell_type,"_sce.RDS"))
  
  ############## Get rough guess of non-DE genes #################
  
  ## Pseudobulk data 
  pseudobulk <- sumCountsAcrossCells(dat_fil, ids = 
                                       DataFrame(sampleID=dat_fil$sampleID,
                                                 subjectID=dat_fil$subjectID,
                                                 time=dat_fil$time)
  )
  
  ## Run model with subject and time as fixed effects
  design_mat=model.matrix(~time+subjectID, data=colData(pseudobulk))
  dge <- edgeR::DGEList(counts = assay(pseudobulk))
  dge<-edgeR::calcNormFactors(dge)
  dge <- edgeR::estimateDisp(dge, design_mat)
  
  ## Hypothesis testing for time effect
  fit_edgeR <- edgeR::glmFit(dge, design_mat)
  my_coef=paste0("time", levels(factor(pseudobulk$time))[2])
  res_tab<-edgeR::glmLRT(fit_edgeR, coef=my_coef)$table
  res_tab$padj=p.adjust(res_tab$PValue, method = "BH")
  
  ##Call genes non-DE if they have padj>0.05 and abs logFC <0.1
  non_de_genes=rownames(res_tab)[res_tab$padj>.05 &abs(res_tab$logFC)<.1]
  
  saveRDS(non_de_genes, paste0("data_processed/non_de_genes/khoo_", my_cell_type,
                               "_non_de_genes.RDS"))
}
