#######################################
##
## Script name: Save khoo data
##
## Purpose of script: Format and save khoo data for simulation
##
## Author: Elizabeth Wynn
##
## Date Created: 2024-09-23
##
#######################################
##
## Notes:
##   
##
#######################################

## load up the packages we will need:
library(Seurat)

#######################################
## 01 - Load in data
#######################################

############## Seurat object downloaded from geo database (GSE196456) #################
load(paste0(myPath, "khoo_data/khoo_raw_dat.Rdata"))

############## Filter to only children #################
dat=Merged_Seuobj_l1[,Merged_Seuobj_l1$Age_group=="Kid"]

############## Format #################
dat$sample=paste(dat$Patient,dat$Status, sep="_")


dat=CreateSeuratObject(counts=dat@assays$RNA@counts,
                       meta.data=dat[[]])


#######################################
## 02 - Save object for each cell types
#######################################

cell_types=c("B", "CD4 T", "CD8 T", "NK")
for(my_cell_type in cell_types){
  ############## Subset to single cell type #################
  print(my_cell_type)
  dat_fil=subset(dat, Final_anno.l1==my_cell_type)
  
  ############## Remove genes that are 0's for all cells #################
  dat_fil=dat_fil[rowSums(dat_fil@assays$RNA$counts)!=0,]
  
  saveRDS(dat_fil, paste0("khoo_data/khoo_", my_cell_type, 
                          ".RDS"))
  
}