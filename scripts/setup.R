# R version: 4.3.0 or later
# Script to check and install only the required packages for the analysis

# 1. Install BiocManager if not installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# 2. Define package lists
cran_packages <- c(
  "ggplot2", "dplyr", "car", "FSA", "rstatix", "reshape2", "tidyverse",
  "tibble", "flextable", "officer", "ggrepel", "tidyr", "iNEXT",
  "RColorBrewer", "ggpubr"
)

bioc_packages <- c(
  "phyloseq", "microbiome", "rcompanion", "PMCMRplus",
  "multcomp", "metagenomeSeq", "DESeq2"
)

# 3. Install CRAN packages if missing
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

# 4. Install Bioconductor packages if missing
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  }
}

message("✔ All specified packages have been checked/installed.")