#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(RTN)
  library(RTNduals) 
  library(arrow)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("CRITICAL: No working directory path provided to R script.")
}
io_dir <- args[1]

# ===== Load Configuration JSON =====
options_file <- file.path(io_dir, "options.json")
if (!file.exists(options_file)) {
  stop("CRITICAL: Configuration options.json missing from workspace directory.")
}
opts <- read_json(options_file, simplifyVector = TRUE)

n_permutations  <- as.integer(opts$n_permutations)
n_bootstrap     <- as.integer(opts$n_bootstrap)
p_adjust_method <- as.character(opts$p_adjust_method)
p_value_cutoff  <- as.numeric(opts$p_value_cutoff)

# ===== Load Staged Datasets =====
matrix_data <- read_parquet(file.path(io_dir, "matrix.parquet"))
tf_vec      <- read_parquet(file.path(io_dir, "tfs.parquet"))[[1]]
target_vec  <- read_parquet(file.path(io_dir, "targets.parquet"))[[1]]

all_cols <- colnames(matrix_data)
index_col_name <- grep("__index_level_0__", all_cols, value = TRUE)
if (length(index_col_name) == 0) index_col_name <- all_cols[1] 

gene_names <- matrix_data[[index_col_name]]
numeric_df <- matrix_data[, !(colnames(matrix_data) %in% index_col_name)]
r_matrix <- as.matrix(numeric_df)
rownames(r_matrix) <- gene_names

valid_targets <- target_vec[target_vec %in% rownames(r_matrix)]

# ===== Run RTN & RTNduals Pipeline =====
cat("Executing network inference pipeline...\n")
tni <- tni.constructor(expData = r_matrix, regulatoryElements = tf_vec)
tni <- tni.permutation(tni, nPermutations = n_permutations, pAdjustMethod = p_adjust_method, globalAdjustment = FALSE)
tni <- tni.bootstrap(tni, nBootstrap = n_bootstrap)
tni <- tni.dpi.filter(tni)

mbr <- tni2mbrPreprocess(tni)
mbr <- mbrAssociation(mbr, pAdjustMethod = p_adjust_method, pValueCutoff = p_value_cutoff)

# ===== 1. Save Complete R State =====
cat("Saving native R data objects...\n")
saveRDS(tni, file = file.path(io_dir, "tni_output.rds"))
saveRDS(mbr, file = file.path(io_dir, "mbr_output.rds"))

# ===== 2. Extract Using Official RTN API Layout =====
cat("Extracting data for pipeline export...\n")
dual_regulons <- mbrGet(mbr, what = "dualRegulons")

grn_df <- data.frame(Regulator1 = character(), Regulator2 = character(), Target = character(), stringsAsFactors = FALSE)

if (!is.null(dual_regulons) && (is.matrix(dual_regulons) || is.data.frame(dual_regulons))) {
  nr <- nrow(dual_regulons)
  if (!is.na(nr) && nr > 0) {
    
    # Retrieve the post-DPI filtered regulons mapping
    regulons <- tni.get(tni, what = "regulons.and.mode")
    
    output_list <- apply(dual_regulons, 1, function(row) {
      tf1 <- row[["Regulator1"]]
      tf2 <- row[["Regulator2"]]
      
      # Extract verified targets for both TFs from the processed regulon list
      targets_tf1 <- names(regulons[[tf1]])
      targets_tf2 <- names(regulons[[tf2]])
      
      # Find valid intersected targets
      shared_tgts <- intersect(targets_tf1, targets_tf2)
      shared_tgts <- intersect(shared_tgts, valid_targets)
      
      if (length(shared_tgts) > 0) {
        return(data.frame(Regulator1 = tf1, Regulator2 = tf2, Target = shared_tgts, stringsAsFactors = FALSE))
      }
      return(NULL)
    }, simplify = FALSE)
    
    output_list <- Filter(Negate(is.null), output_list)
    if (length(output_list) > 0) {
      grn_df <- do.call(rbind, output_list)
    }
  }
}

cat("Total pairwise-target elements mapped:", nrow(grn_df), "\n")
write_parquet(grn_df, file.path(io_dir, "rtnduals_output.parquet"))
cat("--- R Worker Finished ---\n")