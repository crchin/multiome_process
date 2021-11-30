# installation in signac conda environment

# install in conda environment
# conda install -c conda-forge r-hdf5r
# conda install -c bioconda macs2

# Install if needed
install.packages(c("Seurat", "data.table"))
# BiocManager::install("ggbio")
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}

install.packages("tidyverse")

# Install these before installing Signac
BiocManager::install(c("GenomeInfoDb", "GenomicRanges", "IRanges",
                       "Rsamtools", "S4Vectors", "BiocGenerics", "Biostrings", "ggbio", "biovizBase", "AnnotationFilter"))
install.packages("Signac")

#
# BiocManager::install("EnsDb.Mmusculus.v79")
# BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")
# BiocManager::install("SeuratWrappers")
# BiocManager::install(c("EnsDb.Mmusculus.v79", "BSgenome.Mmusculus.UCSC.mm10"))
# BiocManager::install(c("EnsDb.Hsapiens.v86", "BSgenome.Hsapiens.UCSC.hg38", "SeuratWrappers"))

BiocManager::install(
  c(
    "GenomeInfoDb",
    "GenomicRanges",
    "IRanges",
    "Rsamtools",
    "S4Vectors",
    "BiocGenerics",
    "Biostrings",
    "ggbio",
    "biovizBase",
    "AnnotationFilter",
    "EnsDb.Mmusculus.v79",
    "BSgenome.Mmusculus.UCSC.mm10",
    "EnsDb.Hsapiens.v86",
    "BSgenome.Hsapiens.UCSC.hg38"
  )
)


# Make sure hdf5r installed through conda first
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
# remotes::install_github("SeuratDisk")
remotes::install_github("mojaveazure/seurat-disk")
install.packages("R.utils")
remotes::install_github('satijalab/seurat-wrappers')

# For velocity
BiocManager::install("pcaMethods")
devtools::install_github("velocyto-team/velocyto.R")


# usethis::edit_r_environ()
