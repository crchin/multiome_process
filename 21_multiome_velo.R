#! /bin/R
# by Christopher Chin
# 5/2021
# 21 multiome create multiple RNA objects with splicing data from multiome
# Runsheet based

# Packages to install
# Run ahead of time to install packages
# install.packages("Seurat")
# install.packages("data.table")
# BiocManager::install("SeuratWrappers")
# if (!requireNamespace("remotes", quietly = TRUE)) {
#   install.packages("remotes")
# }
# remotes::install_github("SeuratDisk")
# remotes::install_github('satijalab/seurat-wrappers')


# Load data table to read in runsheet (check first)
library(data.table)

# Get runsheet path from arguments
args <- commandArgs(trailingOnly = T)
runsheet.path <- args[1]
# Species must be human or mouse
# species <- args[2]
# # number of cores to use
# if(is.na(args[3])){
#   cores <- 12
# }else{cores <- as.numeric(args[3])}

# Load runsheet and check format
print("Loading runsheet")
runsheet <- fread(runsheet.path, header = T)
# Match runsheet column names
# Path is to cellranger arc outs folder
if (all(c("Sample", "sobj_path") %in% colnames(runsheet)) == T){
  print("Sample and sobj_path column names not found")
  stop()
}else{print("Runsheet accepted!")}

samples <- runsheet$Sample

# Load remaining libraries
library(future)
# plan("multiprocess", workers = cores)
# print(paste("Using", cores, "cores"))

library(Seurat)
library(SeuratWrappers)
library(SeuratDisk)
library(Signac)
library(GenomicRanges)
# Running into memory issues
# library(doMC)
library(velocyto.R)
# registerDoMC(cores)

# Seed set for future
set.seed(42)

options(future.globals.maxSize= 53687091200000000)

# Create directory for intermediate saves
if (!dir.exists("work")){dir.create("work")}
work <- "work"

# Load species based on arg[2]
# if (species %in% c("Human", "human")){
#   library(EnsDb.Hsapiens.v86)
#   library(BSgenome.Hsapiens.UCSC.hg38)
#   annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
#   seqlevelsStyle(annotation) <- "UCSC"
#   genome(annotation) <- "hg38"
#   ref.blacklist <- blacklist_hg38_unified
#   ref.genome <- BSgenome.Hsapiens.UCSC.hg38
# }else if (species %in% c("Mouse","mouse")){
#   library(EnsDb.Mmusculus.v79)
#   library(BSgenome.Mmusculus.UCSC.mm10)
#   annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
#   seqlevelsStyle(annotation) <- "UCSC"
#   genome(annotation) <- "mm10"
#   ref.blacklist <- blacklist_mm10
#   ref.genome <- BSgenome.Mmusculus.UCSC.mm10
# }else {
#   print("species not found")
#   stop()
# }

# Create seurat objects from RNAprint("Creating seurat objects with RNA data")

for (i in samples){
  temp.obj.path <- file.path(work, paste0(i,".sobj.Rdata"))
  temp.outpath <- runsheet$path[runsheet$sample == i]
  
  counts <- Read10X_h5(file.path(temp.outpath, "filtered_feature_bc_matrix.h5"))
  print(paste("Creating",i,"seurat object"))
  temp.sobj <- CreateSeuratObject(
    counts = counts$`Gene Expression`,
    assay = "RNA", project = i
  )
  # temp.sobj <- NormalizeData(temp.sobj, verbose = FALSE)
  # temp.sobj <- FindVariableFeatures(temp.sobj)
  # temp.sobj <- ScaleData(temp.sobj)
  # temp.sobj <- RunPCA(temp.sobj, verbose = T)
  # temp.sobj <- SCTransform(temp.sobj)
  assign(paste0(i,".sobj"), temp.sobj)
  save(list = paste0(i,".sobj"), file = temp.obj.path)
  
  temp.sobj <- RenameCells(temp.sobj, new.names = gsub("-1","",colnames(temp.sobj)))
  
  parent.path <- gsub("outs","",temp.outpath)
  temp.velo.path <- gsub("outs","velocyto",temp.outpath)
  temp.loom <- read.loom.matrices(file.path(temp.velo.path,dir(temp.velo.path)))
  
  for (j in c("spliced","unspliced","ambiguous")){
    colnames(temp.loom[[j]]) <- gsub("x","",gsub(".*:","",colnames(temp.loom[[j]])))
  }
  
  for (j in c("spliced","unspliced","ambiguous")){
    temp <- temp.loom[[j]][,colnames(temp.loom$spliced) %in% colnames(temp.sobj)]
    temp.sobj[[j]] <- CreateAssayObject(temp)
  }
  
  assign(paste0(i,".sobj"), temp.sobj)
  
  print(paste("Saving", paste0(i,".sobj"), "to", file.path(parent.path, paste0(i,"seurat.velocyto.Rdata"))))
  save(list = paste0(i,".sobj"), file = file.path(parent.path, paste0(i,"seurat.velocyto.Rdata")))
  
}

print("All objects created!")



