#! /bin/R
# by Christopher Chin
# 5/2021
# Multiome combine by runsheet
# Combine RNA and combined ATAC peaks

# 08 multiome combine peaks and RNA

# Load and combine single cell data from runsheet
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
  # seurat.velo.path <- file.path(gsub("outs","",temp.outpath),paste0(i,"seurat.velocyto.Rdata"))
  # if(file.exists(seurat.velo.path)){
  #   print(paste("Loading seurat velo object from", seurat.velo.path))
  #   load(seurat.velo.path)
  #   temp.sobj <- get(paste0(i,".sobj"))
  #   counts <- Read10X_h5(file.path(temp.outpath, "filtered_feature_bc_matrix.h5"))
  # }else 
  if(!file.exists(temp.obj.path)){
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

# Merge objects and save
print("Merging multiome objects:")
msobj.list <- paste0(samples, ".msobj")
print(msobj.list)

for (i in msobj.list){
  temp.sobj <- get(i)
  Idents(temp.sobj) <- temp.sobj$orig.ident
  assign(i, temp.sobj)
}

if (length(samples) > 2){
  combined.msobj <- merge(x = get(paste0(samples[1],".msobj")), y = mget(paste0(samples[2:length(samples)], ".msobj")), add.cell.ids = samples)
}else{
  combined.msobj <- merge(x = get(paste0(samples[1],".msobj")), y = get(paste0(samples[2], ".msobj")), add.cell.ids = samples)
}

# Process RNA Data
DefaultAssay(combined.msobj) <- "RNA"
combined.msobj <- NormalizeData(combined.msobj, verbose = FALSE)
combined.msobj <- FindVariableFeatures(combined.msobj)
combined.msobj <- ScaleData(combined.msobj, verbose = T)
combined.msobj <- RunPCA(combined.msobj, npcs = 50, verbose = T)

# Process data
print("Processing combined ATAC data")
DefaultAssay(combined.msobj) <- "peaks"
combined.msobj <- RunTFIDF(combined.msobj)
combined.msobj <- FindTopFeatures(combined.msobj, min.cutoff = "q50")
combined.msobj <- RunSVD(combined.msobj)

print ("Building joint neighbors based on both RNA and ATAC")
combined.msobj <- FindMultiModalNeighbors(
  object = combined.msobj,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:50, 2:40),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)
print("Running UMAP on RNA and ATAC data")
combined.msobj <- RunUMAP(
  object = combined.msobj,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE
)

# link peaks to genes
print("Linking peaks to genes")
DefaultAssay(combined.msobj) <- "peaks"

# first compute the GC content for each peak
combined.msobj <- RegionStats(combined.msobj, genome = ref.genome)

# link peaks to all genes
combined.msobj <- LinkPeaks(
  object = combined.msobj,
  peak.assay = "peaks",
  expression.assay = "RNA"
)

print("Saving Combined multiome object to")
print(file.path(getwd(), "work", "integrated.multiome.sobj.Rdata"))

save(combined.msobj, file = "work/integrated.multiome.sobj.Rdata")
saveRDS(combined.msobj, file = "work/integrated.multiome.sobj.rds")
print("Combined multiome signac object saved!")
