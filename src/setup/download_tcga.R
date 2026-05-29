#!/usr/bin/env Rscript
library(TCGAbiolinks)
library(SummarizedExperiment)

# Accept command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Please provide a project name (e.g., BRCA)", call. = FALSE)
}

project_name <- args[1]
full_project_id <- paste0("TCGA-", project_name)

# Define the base data directory (mounted to the host Mac)
base_dir <- "/data"
output_dir <- file.path(base_dir, project_name)

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

message(paste("Processing project:", full_project_id))

# 1. Query Tumor (added directory parameter)
query_tumor <- GDCquery(
    project = full_project_id,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    sample.type = "Primary Tumor"
)

# 2. Query Normal (added directory parameter)
query_normal <- GDCquery(
    project = full_project_id,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    sample.type = "Solid Tissue Normal"
)

# 3. Download & Prepare (explicitly telling it to download to /data)
GDCdownload(query_tumor, directory = base_dir)
GDCdownload(query_normal, directory = base_dir)

data_tumor <- GDCprepare(query_tumor, directory = base_dir)
data_normal <- GDCprepare(query_normal, directory = base_dir)

matrix_tumor <- assay(data_tumor)
matrix_normal <- assay(data_normal)

# 4. Save to the dynamic folder
write.csv(matrix_tumor, file = file.path(output_dir, paste0(project_name, "_diseased_expression.csv")), row.names = TRUE)
write.csv(matrix_normal, file = file.path(output_dir, paste0(project_name, "_healthy_expression.csv")), row.names = TRUE)

# Print final confirmation and location
message("\n=======================================================")
message("Download complete! Files saved successfully.")
message(paste("Container location: ", output_dir))
message("Your local Mac location: /Users/svj/projects/bio-datasets/", project_name)
message("=======================================================")