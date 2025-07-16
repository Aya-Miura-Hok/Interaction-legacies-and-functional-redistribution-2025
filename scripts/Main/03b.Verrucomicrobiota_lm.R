# Load required libraries
library(phyloseq)
library(microbiome)
library(tidyverse)
library(ggpubr)
library(dplyr)

# Import data
sample_sheet <- read.csv("https://raw.githubusercontent.com/Aya-Miura-Hok/Interaction-legacies-and-functional-redistribution-2025/refs/heads/main/data/bacteria/sample_sheet.csv", row.names = 1)
asv_sheet <- read.csv("https://raw.githubusercontent.com/Aya-Miura-Hok/Interaction-legacies-and-functional-redistribution-2025/refs/heads/main/data/bacteria/otu_table.csv", row.names = 1)
tax_sheet <- read.csv("https://raw.githubusercontent.com/Aya-Miura-Hok/Interaction-legacies-and-functional-redistribution-2025/refs/heads/main/data/bacteria/taxonomy_table.csv", row.names = 1)

# Create phyloseq object
ps_all <- phyloseq(
  otu_table(asv_sheet, taxa_are_rows = FALSE),
  sample_data(sample_sheet),
  tax_table(as.matrix(tax_sheet))
)

# Filter bacterial data and subset treatment samples
ps<-subset_taxa(ps_all, Kingdom == "d__Bacteria")
exp<-subset_samples(ps, Type %in% c("alder_l", "alder_s", "algae_l", "algae_s", "feces_l", "feces_s"))

# Subset samples by treatment
fpomsize  <- subset_samples(exp, Treatment %in% c("C", "I", "M"))
fpomsize1 <- subset_samples(exp, Treatment %in% c("C", "I"))
fpomsize2 <- subset_samples(exp, Treatment %in% c("C", "M"))

# variables for analysis
analyze_fpom <- function(fpom_ps, label) {
  pseq <- fpom_ps %>%
    aggregate_rare(level = "Phylum", detection = 5, prevalence = 0.4) %>%
    microbiome::transform("compositional") %>%
    aggregate_taxa(level = "Phylum")
  
  Verru <- subset_taxa(pseq, Phylum == "Verrucomicrobiota")
  
  df_env <- as(sample_data(Verru), "data.frame") %>%
    rownames_to_column(var = "Sample") %>%
    dplyr::select(Sample, d13C, HIX, BIX)
  
  verru_df <- psmelt(Verru) %>%
    group_by(Sample) %>%
    summarise(Verrucomicrobiota = sum(Abundance)) %>%
    ungroup()
  
  df_lm <- df_env %>%
    left_join(verru_df, by = "Sample") %>%
    drop_na()

  # Linear models and model summaries
  cat(paste0("\n--- ", label, ": d13C ~ Verrucomicrobiota ---\n"))
  print(summary(lm(d13C ~ Verrucomicrobiota, data = df_lm)))
  
  cat(paste0("\n--- ", label, ": HIX ~ Verrucomicrobiota ---\n"))
  print(summary(lm(HIX ~ Verrucomicrobiota, data = df_lm)))
  
  cat(paste0("\n--- ", label, ": BIX ~ Verrucomicrobiota ---\n"))
  print(summary(lm(BIX ~ Verrucomicrobiota, data = df_lm)))
  
  # Visualization
   df_long <- df_lm %>%
    pivot_longer(cols = c(d13C, HIX, BIX), names_to = "Variable", values_to = "Value") %>%
    mutate(Variable = factor(Variable, levels = c("d13C", "HIX", "BIX")))
  
  plot <- ggplot(df_long, aes(x = Verrucomicrobiota, y = Value)) +
    geom_point(alpha = 0.8, size = 2.5) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.9) +
    facet_wrap(~ Variable, scales = "free_y") +
    stat_cor(
      aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
      method = "pearson",
      label.x.npc = "left",
      label.y.npc = "top",
      size = 6
    ) +
    theme_minimal(base_size = 16) +
    labs(
      x = "Verrucomicrobiota relative abundance",
      y = "Response variable",
      title = paste("Linear relationships in", label)
    )
  
  print(plot)
}

# fpom1 (C, I)
fpom1 <- subset_samples(exp, Treatment %in% c("C", "I"))
analyze_fpom(fpom1, "fpom1 (C vs I)")

# fpom2 (C, M)
fpom2 <- subset_samples(exp, Treatment %in% c("C", "M"))
analyze_fpom(fpom2, "fpom2 (C vs M)")
