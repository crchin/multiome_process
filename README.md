# multiome_process
Process and integrate 10x Multiome Data

The following scripts should be used in this order to process the multiome data:

# Create multiome object on cluster from cellranger output path individually (optional)
[01_multiome_signac.R](https://github.com/crchin/multiome_process/blob/main/01_multiome_signac.R)

Input arguments:

argument 1: cellranger out path

argument 2: sample name

argument 3: species, human or mouse


Example Usage:

Rscript ~/scripts/01_multiome_signac.R /athena/melnicklab/users/scratch/chc2077/ls_multiome/ls_multiome/outs ls_multiome mouse

Object will be saved into the cellranger-arc output folder

# Create a csv runsheet that has a unique name for each sample, and a path to cellranger-arc "out" folder
csv runsheet header must be sample path

Example runsheet included:

[20210515_ls_am_multiome_runsheet.csv](https://github.com/crchin/multiome_process/blob/main/20210515_ls_am_multiome_runsheet.csv)

# Load and combine single cell data from runsheet, save as RDS file

This script will create multiome objects if [01_multiome_signac.R](https://github.com/crchin/multiome_process/blob/main/01_multiome_signac.R) was not run

[08_combined_multiome_2.R](https://github.com/crchin/multiome_process/blob/main/08_combined_multiome_2.R)

Input arguments:

argument 1: path to runsheet

argument 2: species to use (Human or Mouse)

argument 3 (optional): number of cores to use

Example Usage:

Rscript ~/.scripts/08_combine_multiome.R ./ls_msobj_combine_runsheet.csv mouse 32

A work directory will be created within the directory that this script is run from. The combined object will be saved as "combined.multiome.sobj.rds" within this work folder

# Use harmony to correct batch effect in both ATAC and RNA assays
[12_multiome_harmony_rds.R](https://github.com/crchin/multiome_process/blob/main/12_multiome_harmony_rds.R)

# Label multiome data by RNA using a reference single cell object
[16_label_singlecell_3.R](https://github.com/crchin/multiome_process/blob/main/16_label_singlecell_3.R)

 
