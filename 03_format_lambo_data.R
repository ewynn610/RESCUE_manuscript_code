#######################################
##
## Script name: Save Lambo data
##
## Purpose of script: Format and save Lambo data for simulation
##
## Author: Elizabeth Wynn
##
#######################################

## load up the packages we will need:

library(Seurat)
library(edgeR)
library(scuttle)
library(Matrix)


#######################################
## 01 -  Read in subject data 
#######################################

# Create a function to load each subject
load_sample <- function(matrix_file) {
  # Derive base name (e.g., GSM7494257_AML16_DX)
  base <- sub("_processed_matrix\\.mtx\\.gz$", "", basename(matrix_file))
  data_dir <- dirname(matrix_file)
  
  # File paths
  matrix_path <- file.path(data_dir, paste0(base, "_processed_matrix.mtx.gz"))
  barcodes_path <- file.path(data_dir, paste0(base, "_processed_barcodes.tsv.gz"))
  genes_path <- file.path(data_dir, paste0(base, "_processed_genes.tsv.gz"))  # or _features
  metadata_path <- file.path(data_dir, paste0(base, "_processed_metadata.tsv.gz"))
  
  # Read matrix and features
  mat <- readMM(matrix_path)
  genes <- read.delim(genes_path, header = FALSE)
  barcodes <- read.delim(barcodes_path, header = FALSE)
  
  # Fix names
  rownames(mat) <- make.unique(genes$V1)  # or genes$V1 depending on format
  colnames(mat) <- barcodes$V1
  
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(counts = mat, project = base)
  
  # Read and add metadata
  metadata <- read.delim(metadata_path)
  metadata <- metadata[match(colnames(seurat_obj), metadata$Cell_Barcode), ]
  rownames(metadata) <- metadata$barcode
  metadata$Cell_Barcode <- NULL 
  
  # Add metadata
  seurat_obj <- AddMetaData(seurat_obj, metadata)
  
  return(seurat_obj)
}


## Define the path to the sample files
## Files include processed data files from GSE235063
sample_files <- list.files(path = paste0("data_raw/lambo_raw/"), pattern = "processed_matrix.mtx.gz", full.names = TRUE)

## Apply to all samples
seurat_list <- lapply(sample_files, load_sample)

#######################################
## 02 - Combine and format data
#######################################


## Merge all into one Seurat object
combined <- merge(seurat_list[[1]], y = seurat_list[-1],
                  add.cell.ids = sapply(seurat_list, function(x) x@project.name))
combined <- JoinLayers(combined)

## Edit labels for consistency
combined$sampleID <- combined$GEO_ID
combined$subjectID <- combined$Lambo_et_al_ID
combined$time <- combined$Patient_Sample

## Remove any samples without all three timepoints
which_rm<-names(which(colSums(table(combined$time, 
                                    combined$subjectID)==0)!=0))

combined=combined[,!combined$subjectID %in% which_rm]


#######################################
## 03 - Save object for each cell types
#######################################

for(x in unique(combined$Classified_Celltype)){
  print(x)
  
  ## Filter data
  combined_fil=combined[,combined$Classified_Celltype==x]
  
  ## Filter to samples with at least 25 cells at all three timepoints
  which_keep<-names(which(table(combined_fil$GEO_ID)>=25))
  combined_fil=combined_fil[,combined_fil$GEO_ID %in% which_keep]
  
  which_keep<-names(which(colSums(table(combined_fil$time, 
                                        combined_fil$subjectID)!=0)==3))
  
  ## Only continue if there are at least 5 subjects left
  if(length(which_keep)>=5){
    combined_fil<-combined_fil[,combined_fil$subjectID %in% which_keep]
    
    ##Filter out genes that are 0 across all cells
    combined_fil=combined_fil[rowSums(combined_fil@assays$RNA$counts)!=0,]
    
    ## Convert to sce
    combined_fil=as.SingleCellExperiment(combined_fil)
    
    ## Save sce object
    saveRDS(combined_fil,paste0("data_processed/sce_objs/lambo_", x, "_sce.RDS"))
    
    
    ############## Get rough guess of non-DE genes #################
    
    ## Pseudobulk data 
    pseudobulk <- sumCountsAcrossCells(combined_fil, ids = 
                                         DataFrame(sampleID=combined_fil$sampleID,
                                                   subjectID=combined_fil$subjectID,
                                                   time=combined_fil$time)
    )
    
    ## Run model with subject and time as fixed effects
    design_mat=model.matrix(~time+subjectID, data=colData(pseudobulk))
    dge <- edgeR::DGEList(counts = assay(pseudobulk))
    dge<-edgeR::calcNormFactors(dge)
    dge <- edgeR::estimateDisp(dge, design_mat)
    fit_edgeR <- edgeR::glmFit(dge, design_mat)
    
    ## Hypothesis testing for pairwise time comparisons
    contrasts<-list(relapse_v_diagnosis=c(0, 1, rep(0, ncol(fit_edgeR$coefficients)-2)),
                    remission_v_diagnosis=c(0,0, 1, rep(0, ncol(fit_edgeR$coefficients)-3)),
                    relapse_v_remission=c(0, 1, -1, rep(0, ncol(fit_edgeR$coefficients)-3)))
    
    ## For each test, keep "non-de" genes with padj>.05 and logFC<.1
    genes_keep=lapply(contrasts, function(x){
      res_tab<-edgeR::glmLRT(fit_edgeR, contrast = x)$table
      res_tab$padj=p.adjust(res_tab$PValue, method = "BH")
      res_tab
      genes_keep=rownames(res_tab)[res_tab$padj>.05&abs(res_tab$logFC)<.1]
    })
    
    ## Keep genes that are non-de for all tests
    non_de_genes=intersect(genes_keep[[1]], intersect(genes_keep[[2]], genes_keep[[3]]))
    
    saveRDS(non_de_genes, paste0("data_processed/non_de_genes/lambo_", x,
                                 "_non_de_genes.RDS"))
  }
}
