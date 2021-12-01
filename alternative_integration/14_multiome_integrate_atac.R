#! /bin/R
# by Christopher Chin
# 5/2021
# 14 Multiome integration
# Seurat/signac integration of ATAC from existing objects

# Takes combined object

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


# Get runsheet path from arguments
args <- commandArgs(trailingOnly = T)
# Object path
msobj.path <- args[1]
# Name of multiome object
msobj.name <- args[2]
# Name of metadata to correct by
batch.name <- args[3]
# number of cores to use
if(is.na(args[4])){
  cores <- 12
}else{cores <- as.numeric(args[4])}

# Load libraries
library(future)
plan("multiprocess", workers = cores)
print(paste("Using", cores, "cores"))

library(Seurat)
library(SeuratWrappers)
library(SeuratDisk)
library(Signac)
library(GenomicRanges)
# Running into memory issues
library(doMC)
registerDoMC(cores)

# Seed set for future
set.seed(42)

options(future.globals.maxSize= 53687091200000000)

print("Loading data")
load(msobj.path)
temp.msobj <- get(msobj.name)

temp.msobj.list <- SplitObject(temp.msobj, split.by = batch.name)

# find integration anchors

print("Finding integration anchors")
integration.anchors <- FindIntegrationAnchors(
  object.list = temp.msobj.list,
  anchor.features = rownames(temp.msobj),
  reduction = "rlsi",
  dims = 2:30
)

# integrate LSI embeddings
print("Integrating LSI embeddings")
temp.atac.int.msobj <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = temp.msobj[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:30
)

temp.atac.int.msobj[["pca"]] <- temp.msobj[["pca"]]

print("Finding multimodal neighbors")
temp.msobj <- FindMultiModalNeighbors(
  object = temp.atac.int.msobj,
  reduction.list = list("pca", "integrated_lsi"), 
  dims.list = list(1:50, 2:30),
  modality.weight.name = c("RNA.weight","peak.weight"),
  verbose = TRUE
)
print("Running UMAP on multimodal data")
temp.msobj <- RunUMAP(
  object = temp.msobj,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE
)

assign(msobj.name, temp.msobj)
print(paste("Saving final object to",file.path(getwd(), paste0(msobj.name, ".RNA.ATAC.integrated.Rdata"))))
save(list = msobj.name, file = file.path(getwd(), paste0(msobj.name, ".RNA.ATAC.integrated.Rdata")))
