#! /bin/R
# by Christopher Chin
# 4/2021
# Multiome ingegrate by runsheet
# Integrate RNA separately

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

# Seed set for future
set.seed(42)

options(future.globals.maxSize= 53687091200)

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
  temp.blacklist <- blacklist_hg38_unified
  temp.genome <- BSgenome.Hsapiens.UCSC.hg38
}else if (species %in% c("Mouse","mouse")){
  library(EnsDb.Mmusculus.v79)
  library(BSgenome.Mmusculus.UCSC.mm10)
  annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
  seqlevelsStyle(annotation) <- "UCSC"
  genome(annotation) <- "mm10"
  temp.blacklist <- blacklist_mm10
  temp.genome <- BSgenome.Mmusculus.UCSC.mm10
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
    temp.sobj <- FindVariableFeatures(temp.sobj)
    temp.sobj <- ScaleData(temp.sobj)
    temp.sobj <- RunPCA(temp.sobj, verbose = T)
    temp.sobj <- SCTransform(temp.sobj)
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

# Integrate RNA data

if (!file.exists("sobj.list.Rdata")){
  sobj.names <- paste0(samples, ".sobj")
  
  print("Objects found:")
  print(sobj.names)
  
  for (i in samples){
    temp.sobj <- get(i)
    DefaultAssay(temp.sobj) <- "RNA"
    temp.sobj[["percent.mt"]] <- PercentageFeatureSet(object = temp.sobj, pattern = "^mt-")
    assign(i, temp.sobj)
    rm(temp.sobj)
  }
  
  # Clear memory
  rm(list = ls(pattern = "\\.sobj"))
  
  sobj.list <- mget(sobj.names)
  
  print("Processing data")
  
  # Process all data
  sobj.list <- lapply(X = sobj.list, FUN = function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, verbose = FALSE)
  })
  
  print("Finding features")
  sobj.features <- SelectIntegrationFeatures(object.list = sobj.list)
  
  print("Scaling data")
  sobj.list <- lapply(X = sobj.list, FUN = function(x) {
    x <- ScaleData(x, features = sobj.features, verbose = FALSE)
    x <- RunPCA(x, features = sobj.features, verbose = FALSE)
  })
  
  print("Saving processed sobj list to")
  print(file.path(getwd(), "work", "sobj.list.Rdata"))
  save(sobj.list, file = "work/sobj.list.Rdata")
  
  
  print("Saving sobj features to")
  print(file.path(getwd(), "work", "sobj.features.Rdata"))
  save(sobj.features, file = "work/sobj.features.Rdata")
}else{
  print("Previous run data found!")
  load("work/sobj.features.Rdata")
  print("sobj.features loaded from")
  print(file.path(getwd(), "work", "sobj.features.Rdata"))
}

if (!file.exists("work/integration.anchors.Rdata")){
  load("work/sobj.list.Rdata")
  print("sobj.list loaded from")
  print(file.path(getwd(), "work", "sobj.list.Rdata"))
  
  print("Finding anchors")
  sobj.anchors <- FindIntegrationAnchors(object.list = sobj.list, reference = c(1, 2), reduction = "rpca", 
                                         dims = 1:50)
  
  print("Saving anchors to")
  print(file.path(getwd(), "work", "integration.anchors.Rdata"))
  save(sobj.anchors, file = "work/integration.anchors.Rdata")
  
  print("Anchors Saved!")
} else{
  print("Previous run data found!")
  load("work/integration.anchors.Rdata")
  print("Anchors loaded from")
  print(file.path(getwd(), "work", "integration.anchors.Rdata"))
}

# Clear memory
rm(sobj.list)

# Integrate objects

print("Integrating object")

int.sobj <- IntegrateData(anchorset = sobj.anchors, dims = 1:50)

int.sobj <- ScaleData(int.sobj, verbose = T)
int.sobj <- RunPCA(int.sobj, npcs = 30, verbose = T)
int.sobj <- RunUMAP(int.sobj, reduction = "pca", dims = 1:30)
int.sobj <- FindNeighbors(object = int.sobj, dims = 1:10)
int.sobj <- FindClusters(object = int.sobj, resolution = 0.5)

print("Saving integrated RNA object to")
print(file.path(getwd(), "work", "integrated.rna.sobj.Rdata"))

save(int.sobj, file = "work/integrated.rna.sobj.Rdata")
print("Integrated RNA seurat object saved!")

# Combine peaks

peak.file.path <- file.path(work, "combined.peaks.Rdata")
if(!file.exists(peak.file.path)){
  for (i in samples){
    temp.sobj <- get(paste0(i,".msobj"))
    print("Calling peaks with MACS2")
    DefaultAssay(temp.sobj) <- "ATAC"
    temp.peaks <- CallPeaks(temp.sobj, macs2.path = "/home/chc2077/miniconda3/envs/signac/bin/macs2")
    temp.peaks <- keepStandardChromosomes(temp.peaks, pruning.mode = "coarse")
    temp.peaks <- subsetByOverlaps(x = temp.peaks, ranges = temp.blacklist, invert = TRUE)
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
    temp.outpath <- runsheet$path[runsheet$sample == i]
    print(paste("Creating",i, "fragment object"))
    temp.frag <- CreateFragmentObject(path = file.path(temp.outpath, "atac_fragments.tsv.gz"))
    assign(paste0(i,".frag"), temp.frag)
    save(list = paste0(i,".frag"), file = file.path(work, paste0(i,".frag.Rdata")))
  }else{
    print(paste("Loading fragment object from", temp.frag.path))
    load(file.path(work, paste0(i,".frag.Rdata")))
  }
}

# Clear objects from RAM
rm(list = ls(pattern = "temp"))

# Create signac objects from combined peaks

for (i in samples){
  temp.msobj <- get(paste0(i,".msobj"))
  temp.msobj.path <- file.path(work, paste0(i,".cp.sobj.Rdata"))
  if(!file.exists(temp.msobj.path)){
    print("Creating ATAC assay and creating object")
    DefaultAssay(temp.msobj) <- "ATAC"
    temp.frag <- get(paste0(i,".frag"))
    temp.counts <- FeatureMatrix(
      fragments = Fragments(temp.msobj),
      features = combined.peaks,
      cells = colnames(temp.msobj)
    )
    
    temp.assay <- CreateChromatinAssay(temp.counts, fragments = temp.frag)
    temp.asobj <- CreateSeuratObject(temp.assay, assay = "peaks")
    
    # Add nucleosome and TSS info to object
    print("Adding Nucleosome signal and TSS enrichment to object")
    DefaultAssay(temp.asobj) <- "peaks"
    temp.asobj <- NucleosomeSignal(temp.asobj)
    temp.asobj <- TSSEnrichment(temp.asobj)
    
    assign(paste0(i,".asobj"), temp.asobj)
    save(list = paste0(i,".asobj"), file = temp.msobj.path)
  }else{
    print(paste("Loading multiome signac object (combined peaks) from", temp.msobj.path))
    load(temp.msobj.path)
  }
}

# Create merged object
print("Merging objects")

signac.objects <- paste0(samples, ".asobj")

if (length(signac.objects) > 2){
  combined <- merge(
    x = get(signac.objects[1]),
    y = mget(signac.objects[2:length(signac.objects)]),
    add.cell.ids = samples
  )
}else{
  combined <- merge(
    x = get(signac.objects[1]),
    y = get(signac.objects[2]),
    add.cell.ids = samples
  )
}

# Add data to integrated RNA



# Process data
print("Processing combined ATAC data")
DefaultAssay(combined) <- "peaks"
combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = "q50")
combined <- RunSVD(combined)

# Print lsi component /sequencing depth plot

pdf(
  file.path(work,"depthcor.pdf")
)
print(
  DepthCor(combined)
)
dev.off()



