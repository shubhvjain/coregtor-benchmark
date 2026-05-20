#!/usr/bin/env Rscript

# Suppress annoying loading messages so Python terminal stays clean
suppressPackageStartupMessages({
  library(CoRegNet)
  library(arrow)
  library(jsonlite)
})

# ===== Process CLI Arguments =====
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("CRITICAL: No working directory path provided to R script.")
}
io_dir <- args[1]

cat("--- R Worker Process Initialized ---\n")
cat("Working directory assigned:", io_dir, "\n")

# ===== Load Configuration JSON =====
options_file <- file.path(io_dir, "options.json")
if (!file.exists(options_file)) {
  stop("CRITICAL: Configuration options.json missing from workspace directory.")
}
opts <- read_json(options_file, simplifyVector = TRUE)

# Parse parameters from the JSON
threshold         <- as.numeric(opts$threshold)
max_coreg         <- as.numeric(opts$max_coreg)
min_coreg_support <- as.numeric(opts$min_coreg_support)
min_gene_support  <- as.numeric(opts$min_gene_support)

# You can easily query your custom extra options dictionary like this:
# my_custom_val  <- opts$extra_config$some_key

cat(sprintf("Parameters - threshold: %s | maxCoreg: %s | minCoregSupport: %s | minGeneSupport: %s\n", 
            threshold, max_coreg, min_coreg_support, min_gene_support))

# ===== Load Staged Datasets =====
cat("Reading data files using Arrow engine...\n")
matrix_data <- read_parquet(file.path(io_dir, "matrix.parquet"))
tf_vec      <- read_parquet(file.path(io_dir, "tfs.parquet"))[[1]]
target_vec  <- read_parquet(file.path(io_dir, "targets.parquet"))[[1]]

# 1. Identify the Index column created by Pandas
# Since get_data doesn't name the index, pyarrow defaults it to "__index_level_0__"
all_cols <- colnames(matrix_data)
index_col_name <- grep("__index_level_0__", all_cols, value = TRUE)

if (length(index_col_name) == 0) {
  # Fallback: check if it used another identifier name
  index_col_name <- all_cols[1] 
}

# 2. Extract the actual gene names vector safely
gene_names <- matrix_data[[index_col_name]]

# 3. Create a numeric matrix containing ONLY the expression measurements
# Drop the index tracking column so we don't pollute the numbers
numeric_df <- matrix_data[, !(colnames(matrix_data) %in% index_col_name)]
r_matrix <- as.matrix(numeric_df)

# 4. Assign the true gene names to the row names of the matrix
rownames(r_matrix) <- gene_names

# ===== Data Validation Check & Sanity Print =====
cat("--- Sanity Check ---\n")
cat("Matrix dimensions:", nrow(r_matrix), "genes x", ncol(r_matrix), "samples\n")
cat("First 3 gene names in matrix:", head(rownames(r_matrix), 3), "\n")
cat("Total TFs requested:", length(tf_vec), " | Total Targets requested:", length(target_vec), "\n")

matched_tfs <- sum(tf_vec %in% rownames(r_matrix))
cat("TFs successfully matched in matrix row names:", matched_tfs, "/", length(tf_vec), "\n")

if (matched_tfs < 2) {
  stop("CRITICAL FAILURE: At least 2 regulators must match the matrix row names. Check serialization indices.")
}
if (anyNA(r_matrix)) {
  stop("CRITICAL FAILURE: Expression matrix contains invalid NA inputs.")
}
cat("--- Sanity Check Passed Perfectly ---\n")

# ===== Discretize =====
cat("Discretizing expression matrix...\n")
disc_matrix <- discretizeExpressionData(r_matrix, threshold = threshold)

# ===== Run hLICORN =====
cat("Running hLICORN...\n")
t0  <- proc.time()
grn <- hLICORN(
  numericalExpression = r_matrix,
  discreteExpression  = disc_matrix,
  TFlist              = tf_vec,
  GeneList            = target_vec,
  parallel            = "no",
  maxCoreg            = max_coreg,
  minCoregSupport     = min_coreg_support,
  minGeneSupport      = min_gene_support,
  verbose             = TRUE
)
hlicorn_secs <- round((proc.time() - t0)[["elapsed"]], 2)

# ===== Export Outputs back to Workspace =====
grn_df <- coregnetToDataframe(grn)
cat("Interactions found:", nrow(grn_df), "\n")
cat("Writing results back to python container mount...\n")
write_parquet(grn_df, file.path(io_dir, "grn_output.parquet"))

cat("--- R Worker Finished Cleanly ---\n")