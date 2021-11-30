#! /bin/R
# by Christopher Chin
# 5/2021
# Label Seurat object saved as RDS
# Make sure seurat RDS file has RDS extension due to how the output file is named

args <- commandArgs(trailingOnly = T)
# Seurat RDS object path
msobj.path <- args[1]
# Path to reference single cell object saved as RDS
ref.sobj.path <- args[2]
# Name of assay to use
assay.name <- args[3]
# Names of labels to transfer
label.names <- args[4:length(args)]

library(Seurat)
library(SeuratWrappers)
library(SeuratDisk)
library(Signac)
library(GenomicRanges)

print("Loading seurat object")
temp.sobj <- readRDS(msobj.path)

print("Loading reference object")
ref.sobj <- readRDS(ref.sobj.path)

print("Setting default assays")
DefaultAssay(temp.sobj) <- assay.name
DefaultAssay(ref.sobj) <- assay.name

temp.list <- list(temp.sobj,ref.sobj)
print("Finding integration features")
anchor.features <- SelectIntegrationFeatures(object.list = temp.list, nfeatures = 2000)

if (!file.exists(gsub("rds","transfer.anchors.rds", msobj.path))){
  print("Finding transfer anchors")
  transfer.anchors <- FindTransferAnchors(reference = ref.sobj, query = temp.sobj, features = anchor.features,
                                          dims = 1:30)
  
  print(paste("Saving transfer anchors to", gsub("rds","transfer.anchors.rds", msobj.path)))
  saveRDS(transfer.anchors, file = gsub("rds","transfer.anchors.rds", msobj.path))
}else{
  transfer.anchors <- readRDS(gsub("rds","transfer.anchors.rds", msobj.path))
}

print("Predicting labels")
for (i in label.names){
  predictions <- TransferData(anchorset = transfer.anchors, refdata = ref.sobj@meta.data[[i]], 
                              dims = 1:30)
  print("Adding labels to object")
  temp.sobj <- AddMetaData(object = temp.sobj, metadata = predictions)
  temp.sobj[[i]] <- temp.sobj$predicted.id
  
}

print(paste("Saving object to", gsub("rds","labeled.rds", msobj.path, ignore.case = T)))
saveRDS(temp.sobj, file = gsub("rds","labeled.rds", msobj.path, ignore.case = T))



