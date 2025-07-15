# Load required packages
library(tidyverse)  # Data handling
library(phyloseq)   # Phylogenetic data analysis
library(microbiome)
library(ggplot2)
library(dplyr)
library(car)        # Levene's test
library(FSA)        # Dunn’s test
library(rstatix)    # Post hoc test for Kruskal–Wallis
library(tibble)
library(flextable)
library(officer)
theme_set(theme_bw())

# Import data from GitHub
sample_sheet <- read.csv("https://raw.githubusercontent.com/Aya-Miura-Hok/Interaction-legacies-and-functional-redistribution-2025/refs/heads/main/data/bacteria/sample_sheet.csv", row.names = 1)
asv_sheet <- read.csv("https://raw.githubusercontent.com/Aya-Miura-Hok/Interaction-legacies-and-functional-redistribution-2025/refs/heads/main/data/bacteria/otu_table.csv", row.names = 1)
tax_sheet <- read.csv("https://raw.githubusercontent.com/Aya-Miura-Hok/Interaction-legacies-and-functional-redistribution-2025/refs/heads/main/data/bacteria/taxonomy_table.csv", row.names = 1)

# Create phyloseq object
ps_all <- phyloseq(otu_table(asv_sheet, taxa_are_rows = FALSE),
                   sample_data(sample_sheet),
                   tax_table(as.matrix(tax_sheet)))

# Filter bacterial samples and FPOM types
ps <- subset_taxa(ps_all, Kingdom == "d__Bacteria")
exp <- subset_samples(ps, Type %in% c("alder_l", "alder_s", "algae_l", "algae_s", "feces_l", "feces_s"))
fpomsize <- subset_samples(exp, Treatment %in% c("C", "I", "M"))

# Aggregate rare taxa and normalize
pseq <- aggregate_rare(fpomsize, level = "Phylum", detection = 5, prevalence = 0.4)
pseq <- microbiome::transform(pseq, transform = "compositional")
pseq <- aggregate_taxa(pseq, level = "Phylum")
pseq_grouped <- merge_samples(pseq, "Type3")
pseq_grouped <- microbiome::transform(pseq_grouped, transform = "compositional")

# Add metadata
sample_info <- strsplit(sample_names(pseq_grouped), "_")
Type_values <- sapply(sample_info, function(x) paste(x[1:(length(x)-1)], collapse = "_"))
Treatment_values <- sapply(sample_info, function(x) x[length(x)])
sample_data(pseq_grouped)$Type <- Type_values
sample_data(pseq_grouped)$Treatment <- Treatment_values

# Convert to dataframe for modeling
df_taxa <- as.data.frame(psmelt(pseq_grouped)) %>%
  dplyr::select(Sample, Phylum, Abundance) %>%
  distinct() %>%
  pivot_wider(names_from = Phylum, values_from = Abundance, values_fill = list(Abundance = 0))

df_env <- as(sample_data(pseq_grouped), "data.frame") %>%
  rownames_to_column(var = "Sample") %>%
  dplyr::select(Sample, d13C, HIX, BIX)

df_corr <- left_join(df_env, df_taxa, by = "Sample")

# Function to run GLM analysis
run_glm_analysis <- function(df_corr, response_var) {
  # 不要な目的変数を除外
  df_filtered <- df_corr %>%
    dplyr::select(-all_of(setdiff(c("d13C", "HIX", "BIX"), response_var)))
  
  # 数値型の説明変数だけに絞って分散の大きいものを抽出
  vars_to_keep <- df_filtered %>%
    dplyr::select(where(is.numeric), -all_of(response_var)) %>%
    summarise_all(var, na.rm = TRUE) %>%
    unlist() %>%
    sort(decreasing = TRUE) %>%
    head(10) %>%
    names()
  
  # 選択された変数だけを残す
  df_reduced <- df_filtered %>%
    dplyr::select(Sample, all_of(response_var), all_of(vars_to_keep))
  
  # GLMとstepwise選択
  model_glm <- glm(as.formula(paste(response_var, "~ .")), data = df_reduced[, -1], family = gaussian())
  model_step <- step(model_glm, direction = "both")
  summary_step <- summary(model_step)
  
  # 決定係数の計算
  rss <- sum(resid(model_step)^2)
  tss <- sum((df_reduced[[response_var]] - mean(df_reduced[[response_var]]))^2)
  r2 <- 1 - (rss / tss)
  
  # 結果の整形
  df_result <- as.data.frame(summary_step$coefficients)
  colnames(df_result) <- c("Estimate", "Std_Error", "t_value", "p_value")
  df_result$Phylum <- rownames(df_result)
  df_result <- df_result %>%
    mutate(Significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ ""
    )) %>%
    filter(Phylum != "(Intercept)") %>%
    mutate(across(c(Estimate, Std_Error, t_value, p_value), round, 2))

  # 表として整形
  ft <- flextable(df_result) %>%
    bold(part = "header") %>%
    align(align = "center", part = "all") %>%
    border_remove() %>%
    hline_top(part = "header", border = fp_border(color = "black", width = 1)) %>%
    hline_bottom(part = "header", border = fp_border(color = "black", width = 1)) %>%
    hline_bottom(part = "body", border = fp_border(color = "black", width = 1)) %>%
    autofit()
  
  return(list(table = ft, r_squared = r2))
}


# Run for each response variable
d13C_results <- run_glm_analysis(df_corr, "d13C")
HIX_results <- run_glm_analysis(df_corr, "HIX")
BIX_results <- run_glm_analysis(df_corr, "BIX")

# View tables
d13C_results$table
HIX_results$table
BIX_results$table