
#! /bin/R
# by Christopher Chin
# 5/2021
# Run harmony on combined object (loaded from RDS)
# Use RDS file as input

# Get runsheet path from arguments
args <- commandArgs(trailingOnly = T)
# RDS Object path
msobj.path <- args[1]
# name of object
msobj.name <- args[2]
# Name of metadata to correct by
batch.name <- args[3]
# number of cores to use
if(is.na(args[4])){
  cores <- 12
}else{cores <- as.numeric(args[4])}

# Check for output object - exit if it exists
if(file.exists(file.path(getwd(), paste0(msobj.name, ".harmony.integrated.rds")))){
  print(paste("Harmony processed seurat objected found at",file.path(getwd(), paste0(msobj.name, ".harmony.integrated.rds"))))
  print("Exiting")
  quit(save = "no")
}


# Load libraries
library(future)
plan("multiprocess", workers = cores)
print(paste("Using", cores, "cores"))

library(Seurat)
library(SeuratWrappers)
library(SeuratDisk)
# Running into memory issues
library(doMC)
registerDoMC(cores)
library(harmony)

# Create work directory
if(!dir.exists("work")){
  print(paste("Creating work directory at", file.path(getwd(),"work")))
  dir.create("work")
}

# Load data
print("Loading data")
temp.msobj <- readRDS(msobj.path)
print("Processing data")
temp.msobj <- FindVariableFeatures(temp.msobj)
temp.msobj <- ScaleData(temp.msobj, verbose = T)
temp.msobj <- RunPCA(temp.msobj, npcs = 50, verbose = T)
print("Running harmony on RNA")
temp.msobj <- RunHarmony(
  object = temp.msobj,
  group.by.vars = batch.name,
  reduction = 'pca',
  assay.use = 'RNA',
  project.dim = FALSE
)
temp.msobj[["harmony_pca"]] <- temp.msobj[["harmony"]]

print("Processing data based on harmony pca")
temp.msobj <- RunUMAP(temp.msobj, reduction = "harmony_pca", dims = 1:30)
temp.msobj <- FindNeighbors(object = temp.msobj, dims = 1:10)
temp.msobj <- FindClusters(object = temp.msobj, resolution = 0.5)

print(paste("Saving harmony corrected object to",file.path(getwd(), paste0(msobj.name, ".harmony.integrated.Rdata"))))

saveRDS(temp.msobj, file = file.path(getwd(), paste0(msobj.name, ".harmony.integrated.rds")))

