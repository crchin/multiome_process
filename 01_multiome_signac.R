#! /bin/R
# by Christopher Chin
# 4/2021
# Create multiome object on cluster from cellranger output path

library(Seurat)
library(Signac)
library(SeuratWrappers)
library(SeuratDisk)
library(data.table)

options(future.globals.maxSize= 53687091200)

# Args:
# args[1] cellranger out path
# args[2] sample name
# args[3] species, human or mouse
args <- commandArgs(trailingOnly = T)

# Load species and assign annotations

if (args[3] %in% c("Human", "human")){
  library(EnsDb.Hsapiens.v86)
  library(BSgenome.Hsapiens.UCSC.hg38)
  annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  seqlevelsStyle(annotation) <- "UCSC"
  genome(annotation) <- "hg38"
  temp.blacklist <- blacklist_hg38_unified
  temp.genome <- BSgenome.Hsapiens.UCSC.hg38
}else if (args[3] %in% c("Mouse","mouse")){
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

# Create object
print(paste("Loading counts from", file.path(args[1], "filtered_feature_bc_matrix.h5")))
counts <- Read10X_h5(file.path(args[1], "filtered_feature_bc_matrix.h5"))
fragpath <- file.path(args[1], "atac_fragments.tsv.gz")

# create a Seurat object containing the RNA data
print("Creating seurat object with RNA data")
multiome.sobj <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

# create ATAC assay and add it to the object
print("Creating ATAC assay and adding to object")
multiome.sobj[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)

print("multiome object:")
multiome.sobj

# Set default assay to ATAC

DefaultAssay(multiome.sobj) <- "ATAC"

# Add nucleosome and TSS info to object
print("Adding Nucleosome signal and TSS enrichment to object")
multiome.sobj <- NucleosomeSignal(multiome.sobj)
multiome.sobj <- TSSEnrichment(multiome.sobj)

# Call peaks with Macs2
print("Calling peaks with MACS2")
peaks <- CallPeaks(multiome.sobj, macs2.path = "/home/chc2077/miniconda3/envs/signac/bin/macs2")

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = temp.blacklist, invert = TRUE)

# quantify counts in each peak
print("Quantifying counts per peak")
macs2_counts <- FeatureMatrix(
  fragments = Fragments(multiome.sobj),
  features = peaks,
  cells = colnames(multiome.sobj)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
print("Adding peak data to object")
multiome.sobj[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragpath,
  annotation = annotation
)


# process RNA data
print ("Processing RNA data")
DefaultAssay(multiome.sobj) <- "RNA"
multiome.sobj <- SCTransform(multiome.sobj)
multiome.sobj <- RunPCA(multiome.sobj)


# process macs2 peak data
print ("Processing peak data")
DefaultAssay(multiome.sobj) <- "peaks"
multiome.sobj <- FindTopFeatures(multiome.sobj, min.cutoff = 5)
multiome.sobj <- RunTFIDF(multiome.sobj)
multiome.sobj <- RunSVD(multiome.sobj)


# build a joint neighbor graph using both assays
print ("Building joint neighbors based on both RNA and ATAC")
multiome.sobj <- FindMultiModalNeighbors(
  object = multiome.sobj,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:50, 2:40),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)

# build a joint UMAP visualization
print("Running UMAP on RNA and ATAC data")
multiome.sobj <- RunUMAP(
  object = multiome.sobj,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE
)

# link peaks to genes
print("Linking peaks to genes")
DefaultAssay(multiome.sobj) <- "peaks"

# first compute the GC content for each peak
multiome.sobj <- RegionStats(multiome.sobj, genome = temp.genome)

# link peaks to all genes
multiome.sobj <- LinkPeaks(
  object = multiome.sobj,
  peak.assay = "peaks",
  expression.assay = "SCT"
)

print("Saving object")
save(multiome.sobj, file = paste0(args[2], ".processed.multiome.sobj.Rdata"))
saveRDS(multiome.sobj, file = paste0(args[2], ".processed.multiome.sobj.RDS"))