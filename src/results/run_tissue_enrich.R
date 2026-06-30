args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop("Usage: Rscript run_tissue_enrich.R <project_root> <input_csv> <output_json>")
}

project_root <- args[1]
input_csv <- args[2]
output_json <- args[3]

suppressPackageStartupMessages(library(renv))
renv::load(project_root)

suppressPackageStartupMessages(library(TissueEnrich))
suppressPackageStartupMessages(library(GSEABase))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(jsonlite))

coreg_file <- read.csv(input_csv, stringsAsFactors = FALSE)

if (nrow(coreg_file) == 0) {
  stop("Input CSV has no rows")
}

results <- vector("list", nrow(coreg_file))

for (i in seq_len(nrow(coreg_file))) {
  row <- coreg_file[i, ]

  genes <- strsplit(row$common_targets, ";")[[1]]
  genes <- trimws(genes)
  genes <- unique(genes[genes != ""])

  if (length(genes) == 0) {
    results[[i]] <- list(
      row_index = i,
      source = row$source,
      support_size = unname(row$support_size),
      note = row$note,
      genes = character(0),
      tissue_results = list(),
      error = "No genes found in common_targets"
    )
    next
  }

  input_gs <- GeneSet(
    geneIds = genes,
    organism = "Homo Sapiens",
    geneIdType = SymbolIdentifier()
  )

  res <- tryCatch(
    TissueEnrich::teEnrichment(
      inputGenes = input_gs,
      rnaSeqDataset = 2
    ),
    error = function(e) e
  )

  if (inherits(res, "error")) {
    results[[i]] <- list(
      row_index = i,
      source = row$source,
      support_size = unname(row$support_size),
      note = row$note,
      genes = genes,
      tissue_results = list(),
      error = conditionMessage(res)
    )
    next
  }

  summary_se <- res[[1]]
  mat <- SummarizedExperiment::assay(summary_se)
  df <- as.data.frame(mat, stringsAsFactors = FALSE)
  df$tissue <- rownames(mat)
  df <- df[, c("tissue", setdiff(names(df), "tissue"))]

  tissue_results <- apply(df, 1, function(x) {
    row_list <- as.list(x)

    for (nm in names(row_list)) {
      if (nm != "tissue" && !is.null(row_list[[nm]]) && !is.na(row_list[[nm]])) {
        suppressWarnings({
          num_val <- as.numeric(row_list[[nm]])
          if (!is.na(num_val)) {
            row_list[[nm]] <- num_val
          }
        })
      }
    }

    row_list
  })

  tissue_results <- unname(as.list(tissue_results))

  results[[i]] <- list(
    row_index = i,
    source = row$source,
    support_size = unname(row$support_size),
    note = row$note,
    genes = genes,
    tissue_results = tissue_results
  )
}

payload <- list(
  input_rows = nrow(coreg_file),
  dataset = "GTEx",
  organism = "Homo Sapiens",
  results = results
)

write_json(payload, output_json, pretty = TRUE, auto_unbox = TRUE, na = "null")