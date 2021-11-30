# multiome_process
Process and integrate 10x Multiome Data

The following scripts should be used in this order to process the multiome data:

# Create multiome object on cluster from cellranger output path
01_multiome_signac.R

# Create a csv runsheet that has a unique name for each sample, and a path to cellrager "out" folder
# csv runsheet header must be sample path
# Example runsheet included
20210515_ls_am_multiome_runsheet.csv

# Load and combine single cell data from runsheet, save as RDS file
08_combined_multiome_2.R

# Use harmony to correct batch effect in both ATAC and RNA assays
12_multiome_harmony_rds.R

# Label multiome data by RNA using a reference single cell object
16_label_singlecell_3.R

