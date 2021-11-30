#! /bin/R
# by Christopher Chin
# 5/2021
# Convert Seurat object saved as RData to RDS

library(Seurat)
library(SeuratWrappers)
library(SeuratDisk)
library(Signac)
library(GenomicRanges)

# Get runsheet path from arguments
args <- commandArgs(trailingOnly = T)
# Path to Rdata
msobj.path <- args[1]
# Name of object
msobj.name <- args[2]
# output path
out.path <- args[3]

load(msobj.path)
temp.object <- get(msobj.name)

saveRDS(temp.object, file = out.path)
