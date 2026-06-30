# setup_renv.R

# 0. Ensure renv is available and activate it
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv", repos = "https://cloud.r-project.org")
}
renv::activate()

# 1. Install ensurer 1.1 from CRAN archive via renv
if (!requireNamespace("ensurer", quietly = TRUE) || packageVersion("ensurer") < "1.1") {
  renv::install("https://cran.r-project.org/src/contrib/Archive/ensurer/ensurer_1.1.tar.gz")
}

# 2. Install TissueEnrich via renv (Bioconductor)
if (!requireNamespace("TissueEnrich", quietly = TRUE)) {
  renv::install("TissueEnrich")
}

if (!requireNamespace("future.apply", quietly = TRUE)) {
  renv::install("future.apply")
}

# 3. Snapshot and verify
renv::snapshot()
library(TissueEnrich)
library(future)
library(future.apply)
cat("Setup complete inside renv! TissueEnrich and dependencies loaded successfully.\n")