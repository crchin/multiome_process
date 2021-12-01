#! /bin/R
# by Christopher Chin
# 5/2021
# 13 Multiome integration
# Seurat/signac integration of both RNA and ATAC
# Didn't work - fragment name issues?
# Need to change order


# Load and integrate single cell data from runsheet
# Included if portions so analysis can be restarted

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
species <- args[2]
# number of cores to use
if(is.na(args[3])){
  cores <- 12
}else{cores <- as.numeric(args[3])}

# Load runsheet and check format
print("Loading runsheet")
runsheet <- fread(runsheet.path, header = T)
# Match runsheet column names
# Path is to cellranger arc outs folder
if (all(colnames(runsheet) == c("sample","path")) == F){
  print("Runsheet colnames do not match")
  stop()
}else{print("Runsheet accepted!")}

samples <- runsheet$sample

# Load remaining libraries
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

# Create directory for intermediate saves
if (!dir.exists("work")){dir.create("work")}
work <- "work"

# Load species based on arg[2]
if (species %in% c("Human", "human")){
  library(EnsDb.Hsapiens.v86)
  library(BSgenome.Hsapiens.UCSC.hg38)
  annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  seqlevelsStyle(annotation) <- "UCSC"
  genome(annotation) <- "hg38"
  ref.blacklist <- blacklist_hg38_unified
  ref.genome <- BSgenome.Hsapiens.UCSC.hg38
}else if (species %in% c("Mouse","mouse")){
  library(EnsDb.Mmusculus.v79)
  library(BSgenome.Mmusculus.UCSC.mm10)
  annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
  seqlevelsStyle(annotation) <- "UCSC"
  genome(annotation) <- "mm10"
  ref.blacklist <- blacklist_mm10
  ref.genome <- BSgenome.Mmusculus.UCSC.mm10
}else {
  print("species not found")
  stop()
}

# Create seurat objects from RNA
print("Creating seurat objects with RNA data")

for (i in samples){
  temp.obj.path <- file.path(work, paste0(i,".sobj.Rdata"))
  temp.outpath <- runsheet$path[runsheet$sample == i]
  if(!file.exists(temp.obj.path)){
    counts <- Read10X_h5(file.path(temp.outpath, "filtered_feature_bc_matrix.h5"))
    print(paste("Creating",i,"seurat object"))
    temp.sobj <- CreateSeuratObject(
      counts = counts$`Gene Expression`,
      assay = "RNA", project = i
    )
    
    temp.sobj <- NormalizeData(temp.sobj, verbose = FALSE)
    temp.sobj <- FindVariableFeatures(temp.sobj)
    temp.sobj <- ScaleData(temp.sobj)
    temp.sobj <- RunPCA(temp.sobj, verbose = T)
    # temp.sobj <- SCTransform(temp.sobj)
    assign(paste0(i,".sobj"), temp.sobj)
    save(list = paste0(i,".sobj"), file = temp.obj.path)
  }else{
    print(paste("Loading seurat object from", temp.obj.path))
    load(temp.obj.path)
    temp.sobj <- get(paste0(i,".sobj"))
  }
  temp.msobj.path <- file.path(work, paste0(i,".indv.peaks.msobj.Rdata"))
  if(!file.exists(temp.msobj.path)){
    print("Creating ATAC assay and adding to object")
    temp.sobj[["ATAC"]] <- CreateChromatinAssay(
      counts = counts$Peaks,
      sep = c(":", "-"),
      fragments = file.path(temp.outpath, "atac_fragments.tsv.gz"),
      annotation = annotation
    )
    assign(paste0(i,".msobj"), temp.sobj)
    save(list = paste0(i,".msobj"), file = temp.msobj.path)
  }else{
    print(paste("Loading multiome seurat object (individual peaks) from", temp.msobj.path))
    load(temp.msobj.path)
  }
}


# Combine peaks

peak.file.path <- file.path(work, "combined.peaks.Rdata")
if(!file.exists(peak.file.path)){
  for (i in samples){
    temp.sobj <- get(paste0(i,".msobj"))
    print("Calling peaks with MACS2")
    DefaultAssay(temp.sobj) <- "ATAC"
    temp.peaks <- CallPeaks(temp.sobj, macs2.path = "/home/chc2077/miniconda3/envs/signac/bin/macs2")
    temp.peaks <- keepStandardChromosomes(temp.peaks, pruning.mode = "coarse")
    temp.peaks <- subsetByOverlaps(x = temp.peaks, ranges = ref.blacklist, invert = TRUE)
    assign(paste0(i,".peaks"), temp.peaks)
  }
  print("Combining peaks")
  all.peaks <- paste0(samples, ".peaks")
  combined.peaks <- get(all.peaks[1])
  for (i in all.peaks[2:length(all.peaks)]){
    combined.peaks <- c(combined.peaks,get(i))
  }
  combined.peaks <- reduce(combined.peaks)
  save(combined.peaks, file = peak.file.path)
}else{
  print(paste("Loading peaks from", peak.file.path))
  load(peak.file.path)
}

# Filter peaks
print("Flitering peaks")
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]

# Create fragment objects
for (i in samples){
  temp.frag.path <- file.path(work, paste0(i,".frag.Rdata"))
  if(!file.exists(temp.frag.path)){
    temp.sobj <- get(paste0(i,".sobj"))
    temp.outpath <- runsheet$path[runsheet$sample == i]
    print(paste("Creating",i, "fragment object"))
    temp.frag <- CreateFragmentObject(path = file.path(temp.outpath, "atac_fragments.tsv.gz"), cells = colnames(temp.sobj))
    assign(paste0(i,".frag"), temp.frag)
    save(list = paste0(i,".frag"), file = file.path(work, paste0(i,".frag.Rdata")))
  }else{
    print(paste("Loading fragment object from", temp.frag.path))
    load(file.path(work, paste0(i,".frag.Rdata")))
  }
}


# Clear objects from RAM
rm(list = ls(pattern = "temp"))

# Add ATAC data for combined peaks
for (i in samples){
  temp.msobj <- get(paste0(i,".msobj"))
  temp.sobj <- get(paste0(i,".sobj"))
  temp.msobj.path <- file.path(work, paste0(i,".combined.peaks.msobj.Rdata"))
  if(!file.exists(temp.msobj.path)){
    print("Creating ATAC assay and adding to object")
    DefaultAssay(temp.msobj) <- "ATAC"
    temp.frag <- get(paste0(i,".frag"))
    temp.counts <- FeatureMatrix(
      fragments = Fragments(temp.msobj),
      features = combined.peaks,
      cells = colnames(temp.msobj)
    )
    temp.sobj[["peaks"]] <- CreateChromatinAssay(
      counts = temp.counts,
      fragments = temp.frag,
      annotation = annotation
    )
    # Add nucleosome and TSS info to object
    print("Adding Nucleosome signal and TSS enrichment to object")
    DefaultAssay(temp.sobj) <- "peaks"
    temp.sobj <- NucleosomeSignal(temp.sobj)
    temp.sobj <- TSSEnrichment(temp.sobj)
    
    temp.sobj <- RenameCells(temp.sobj, new.names = paste(i, colnames(temp.sobj), sep = "_"))
    
    assign(paste0(i,".msobj"), temp.sobj)
    save(list = paste0(i,".msobj"), file = temp.msobj.path)
  }else{
    print(paste("Loading multiome seurat object (combined peaks) from", temp.msobj.path))
    load(temp.msobj.path)
  }
}

# Integrate RNA data
sobj.names <- paste0(samples, ".msobj")

if (!file.exists("work/sobj.list.Rdata")){
  
  print("Objects found:")
  print(sobj.names)
  
  for (i in sobj.names){
    temp.sobj <- get(i)
    DefaultAssay(temp.sobj) <- "RNA"
    Idents(temp.sobj) <- temp.sobj$orig.ident
    temp.sobj[["percent.mt"]] <- PercentageFeatureSet(object = temp.sobj, pattern = "^mt-")
    assign(i, temp.sobj)
    rm(temp.sobj)
  }
  # sobj.list <- mget(sobj.names)
  
  sobj.list <- list()
  for (i in sobj.names){
    sobj.list[[i]] <- get(i)
  }
  
  # Clear memory
  rm(list = ls(pattern = "\\.sobj"))
  
  # print("Processing data")
  # Process all data
  # sobj.list <- lapply(X = sobj.list, FUN = function(x) {
  #   x <- NormalizeData(x, verbose = FALSE)
  #   x <- FindVariableFeatures(x, verbose = FALSE)
  # })
  
  print("Finding features")
  sobj.features <- SelectIntegrationFeatures(object.list = sobj.list)
  
  print("Saving sobj features to")
  print(file.path(getwd(), "work", "sobj.features.Rdata"))
  save(sobj.features, file = "work/sobj.features.Rdata")
  
  print("Scaling data")
  sobj.list <- lapply(X = sobj.list, FUN = function(x) {
    x <- ScaleData(x, features = sobj.features, verbose = FALSE)
    x <- RunPCA(x, features = sobj.features, verbose = FALSE)
  })
  
  print("Saving processed sobj list to")
  print(file.path(getwd(), "work", "sobj.list.Rdata"))
  save(sobj.list, file = "work/sobj.list.Rdata")
  
}else{
  print("Previous run data found!")
  load("work/sobj.features.Rdata")
  print("sobj.features loaded from")
  print(file.path(getwd(), "work", "sobj.features.Rdata"))
}

if (!file.exists("work/integration.anchors.Rdata")){
  if (!exists("sobj.list")){
    print("Loading sobj.list")
    load("work/sobj.list.Rdata")
    print("sobj.list loaded from")
    print(file.path(getwd(), "work", "sobj.list.Rdata"))
  }
  print("Removing individual sobjs")
  rm(list = ls(pattern = "\\.msobj"))
  rm(list = ls(pattern = "\\.sobj"))
  
  for (i in sobj.names){
    Idents(sobj.list[[i]]) <- sobj.list[[i]]$orig.ident
  }
  
  
  print("Finding anchors")
  if (length(sobj.list) > 2){
    print("Integrating by rpca")
    sobj.anchors <- FindIntegrationAnchors(object.list = sobj.list, reference = c(1, 2), reduction = "rpca", 
                                           dims = 1:50)
  }else{
    print("Using conventional integration")
    sobj.anchors <- FindIntegrationAnchors(object.list = sobj.list, normalization.method = "LogNormalize",
                                           anchor.features = sobj.features, verbose = T)
  }
  
  print("Saving anchors to")
  print(file.path(getwd(), "work", "integration.anchors.Rdata"))
  save(sobj.anchors, file = "work/integration.anchors.Rdata")
  
  print("Anchors Saved!")
} else{
  # print("Previous run data found!")
  print("Loading anchors from previous run")
  load("work/integration.anchors.Rdata")
  print("Anchors loaded from")
  print(file.path(getwd(), "work", "integration.anchors.Rdata"))
}

# Clear memory
if(exists("sobj.list")){
  rm(sobj.list)
}

# Integrate objects

if (!file.exists("work/unprocessed.int.sobj.Rdata")){
  print("Integrating object")
  int.sobj <- IntegrateData(anchorset = sobj.anchors, dims = 1:50)
  print("Saving unprocessed integrated object to")
  print(file.path(getwd(), "work", "unprocessed.int.sobj.Rdata"))
  save(int.sobj, file = "work/unprocessed.int.sobj.Rdata")
}else{
  print("Loading integrated object from previous run")
  load("work/unprocessed.int.sobj.Rdata")
}

int.sobj <- ScaleData(int.sobj, verbose = T)
int.sobj <- RunPCA(int.sobj, npcs = 50, verbose = T)

# Process data
print("Processing combined ATAC data")
DefaultAssay(int.sobj) <- "peaks"
int.sobj <- RunTFIDF(int.sobj)
int.sobj <- FindTopFeatures(int.sobj, min.cutoff = "q50")
int.sobj <- RunSVD(int.sobj)

# Print lsi component /sequencing depth plot

pdf(
  file.path(work,"depthcor.pdf")
)
print(
  DepthCor(int.sobj)
)
dev.off()

print ("Building joint neighbors based on both RNA and ATAC")
int.sobj <- FindMultiModalNeighbors(
  object = int.sobj,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:50, 2:40),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)
print("Running UMAP on RNA and ATAC data")
int.sobj <- RunUMAP(
  object = int.sobj,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE
)

# int.sobj <- RunUMAP(int.sobj, reduction = "pca", dims = 1:30)
# int.sobj <- FindNeighbors(object = int.sobj, dims = 1:10)
# int.sobj <- FindClusters(object = int.sobj, resolution = 0.5)
# DefaultAssay(int.sobj) <- "integrated"
# int.sobj <- FindClusters(object = int.sobj, verbose = FALSE, algorithm = 3)

# link peaks to genes
print("Linking peaks to genes")
DefaultAssay(int.sobj) <- "peaks"

# first compute the GC content for each peak
int.sobj <- RegionStats(int.sobj, genome = ref.genome)

# link peaks to all genes
int.sobj <- LinkPeaks(
  object = int.sobj,
  peak.assay = "peaks",
  expression.assay = "RNA"
)


print("Saving RNA integrated multiome object to")
print(file.path(getwd(), "work", "integrated.multiome.sobj.Rdata"))

save(int.sobj, file = "work/integrated.multiome.sobj.Rdata")
print("RNA Integrated multiome signac object saved!")




