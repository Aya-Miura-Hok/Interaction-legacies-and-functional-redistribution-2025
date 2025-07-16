# Load required libraries
library(ggplot2)
library(microbiome)
library(phyloseq)
library(dplyr)
library(ggpubr)
library(tidyr)
library(tibble)
library(tidyverse)
library(RColorBrewer)

# Set default ggplot theme
theme_set(theme_bw())

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
fpomsize <- subset_samples(exp, Treatment %in% c("C", "I", "M"))

# Aggregate rare taxa and normalize by relative abundance
pseq_phylum <- fpomsize %>%
  aggregate_rare(level = "Phylum", detection = 5, prevalence = 0.2) %>%
  microbiome::transform("compositional") %>%
  aggregate_taxa("Phylum") %>%
  merge_samples("Type3") %>%
  microbiome::transform("compositional")

# Focus on Verrucomicrobiota
Verru <- subset_taxa(pseq_phylum, Phylum == "Verrucomicrobiota")

# Extract relative abundance and summarize per sample
df_phylum <- psmelt(Verru) %>%
  group_by(Sample, Phylum) %>%
  summarize(Abundance = sum(Abundance), .groups = "drop")

# Derive 'Treatment' label from the last character of sample name
df_phylum <- df_phylum %>%
  mutate(Treatment = substr(Sample, nchar(Sample), nchar(Sample)))

# Convert 'Treatment' to factor and ensure correct order
df_phylum$Treatment <- factor(df_phylum$Treatment, levels = c("C", "I", "M"))

# --- Statistical Analysis ---

# ANOVA
anova_model <- aov(Abundance ~ Treatment, data = df_phylum)
print(summary(anova_model))

# Normality & homogeneity of variance checks
shapiro.test(residuals(anova_model))  # Normality
car::leveneTest(Abundance ~ Treatment, data = df_phylum)  # Homoscedasticity

# Kruskal-Wallis test as non-parametric alternative
print(kruskal.test(Abundance ~ Treatment, data = df_phylum))

# Tukey's HSD post hoc test
tukey_result <- TukeyHSD(anova_model)
print(tukey_result)

# Pairwise Wilcoxon test with Bonferroni correction
pairwise_result <- pairwise.wilcox.test(df_phylum$Abundance, df_phylum$Treatment, p.adjust.method = "bonferroni")
print(pairwise_result)

# Format Tukey HSD results for annotation
tukey_df <- as.data.frame(tukey_result$Treatment) %>%
  rownames_to_column("Comparison") %>%
  separate(Comparison, into = c("group1", "group2"), sep = "-") %>%
  mutate(
    significance = case_when(
      `p adj` < 0.001 ~ "***",
      `p adj` < 0.01  ~ "**",
      `p adj` < 0.05  ~ "*",
      TRUE ~ ""
    ),
    p_value = `p adj`
  ) %>%
  filter(significance != "")  # Remove non-significant comparisons

# Convert group labels to factors
tukey_df$group1 <- factor(tukey_df$group1, levels = c("C", "I", "M"))
tukey_df$group2 <- factor(tukey_df$group2, levels = c("C", "I", "M"))

# Define y-axis positions to avoid label overlap
tukey_df <- tukey_df %>%
  mutate(y.position = seq(from = max(df_phylum$Abundance) * 1.2, 
                          by = 0.02, length.out = nrow(tukey_df)))

# Final plot: Boxplot with significance annotations
p <- ggplot(df_phylum, aes(x = Treatment, y = Abundance)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +
  stat_compare_means(method = "anova", label.y = max(df_phylum$Abundance) * 1.1) +
  stat_pvalue_manual(tukey_df, label = "significance", size = 5) +
  theme_minimal() +
  labs(
    title = "Verrucomicrobiota Abundance by Treatment",
    x = "Treatment",
    y = "Relative Abundance"
  )

print(p)
