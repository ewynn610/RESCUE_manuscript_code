#######################################
##
## Script name: Save Mould Data
##
## Purpose of script: Format and save Mould data for simulations
##
## Author: Elizabeth Wynn
#######################################

## load up the packages we will need:
library(Seurat)
library(edgeR)
library(scuttle)

#######################################
## 01 - Read in and filter data
#######################################

## Read in Seurat object containing data
dat_full=readRDS("data_raw/mould_raw.RDS")

## Make new seurat object with only useful information
dat_full=CreateSeuratObject(counts=dat_full@assays$RNA@counts,
                            meta.data = data.frame(cell_type=dat_full$overall_celltype,
                                                   sampleID=dat_full$sample_number,
                                                   subjectID=dat_full$subject,
                                                   time=dat_full$time))

## Only keep subjects with with days 4 and 5 follow-up
subj_ids_keep=c(1714, 3156, 4683, 6386, 28015)
dat_fil=subset(dat_full, subjectID%in%subj_ids_keep)

## Define time as baseline or follow-up
dat_fil$time=ifelse(dat_fil$time==0, "baseline", "follow-up")

#######################################
## 02 - Save data for each celltype
#######################################

## Get all cell tiypes
cell_types=unique(dat_fil$cell_type)

##Loop through celltypes
for(my_cell_type in cell_types){
  print(my_cell_type)
  
  ## Subset down to single cell type 
  dat=subset(dat_fil, cell_type==my_cell_type)
  
  ## If < 25 cells in any sample skip
  ## We want at least 25 cells for both timepoints across >=5 subjects
  ## We only have 5 subjs
  if(any(table(dat$sampleID)<25)){
    next
  }

  ## Remove genes that have 0's across all cells
  genes_keep=which(rowSums(dat@assays$RNA$counts)!=0)
  dat=dat[genes_keep,]
  
  ## Convert to SCE and save
  dat<-as.SingleCellExperiment(dat)
  saveRDS(dat, paste0("data_processed/sce_objs/mould_", my_cell_type, "_sce.RDS"))
  
  ############## Get rough guess of non-DE genes #################
  
  ## Pseudobulk data 
  pseudobulk <- sumCountsAcrossCells(dat, ids = 
                                       DataFrame(sampleID=dat$sampleID,
                                                 subjectID=dat$subjectID,
                                                 time=dat$time)
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
  
  saveRDS(non_de_genes, paste0("data_processed/non_de_genes/mould_", my_cell_type, "_non_de_genes.RDS"))
}

