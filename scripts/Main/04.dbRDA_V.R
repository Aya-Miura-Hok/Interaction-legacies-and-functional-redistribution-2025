# Load required libraries
library(tidyverse)
library(phyloseq)
library(ggrepel)
library(ggplot2)
library(microbiome)
library(vegan)
library(dplyr)
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
ps <- subset_taxa(ps_all, Kingdom == "d__Bacteria")
exp <- subset_samples(ps, Type %in% c("alder_l", "alder_s", "algae_l", "algae_s", "feces_l", "feces_s"))
fpomsize <- subset_samples(exp, Treatment %in% c("C", "I", "M"))

# Aggregate rare phyla and apply compositional transformation
pseq <- aggregate_rare(fpomsize, level = "Phylum", detection = 5, prevalence = 0.2)
pseq <- microbiome::transform(pseq, transform = "compositional")
pseq <- aggregate_taxa(pseq, level = "Phylum")
pseq <- microbiome::transform(pseq, transform = "compositional")

# Prepare OTU and taxonomy tables
otu_table_raw <- t(as.data.frame(otu_table(pseq)))
taxonomy <- as.data.frame(tax_table(pseq))
taxonomy$Sample <- rownames(taxonomy)

otu_table_long <- otu_table_raw %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "Sample") %>%
    left_join(taxonomy, by = "Sample")

# Select key phyla(selected by GLM analysis)
selected_phyla <- c("Proteobacteria", "Verrucomicrobiota", "Planctomycetota", "Desulfobacterota", "Bacteroidota", "Cyanobacteria")
df_phylum_filtered <- dplyr::select(otu_table_long, Sample, all_of(selected_phyla))

# Prepare environmental data
df_env <- as(sample_data(fpomsize), "data.frame") %>%
  rownames_to_column(var = "Sample") %>%
  dplyr::select(Sample, Type, Treatment, d15N, d13C, HIX, BIX, DON)

# Merge OTU and environmental data
df_combined <- inner_join(df_phylum_filtered, df_env, by = "Sample")
df_phylum <- df_combined %>% dplyr::select(Sample, all_of(selected_phyla)) %>% column_to_rownames(var = "Sample")
df_envs <- df_combined %>% dplyr::select(Sample, d13C, HIX, BIX, DON) %>% column_to_rownames(var = "Sample")

# Clean and transform data
df_phylum_fil <- df_phylum %>% filter(complete.cases(.))
df_phylum_fil <- df_phylum_fil[rowSums(df_phylum_fil) > 0, ]
df_phylum_hel <- decostand(df_phylum_fil, method = "hellinger")
bray_dist <- vegdist(df_phylum_hel, method = "bray")
df_envs <- df_envs[rownames(df_phylum_hel), , drop = FALSE]
df_envs_scaled <- scale(df_envs)

# dbRDA analysis
dbrda_model <- capscale(bray_dist ~ ., data = as.data.frame(df_envs_scaled))
summary(dbrda_model)
anova.cca(dbrda_model, permutations = 999)
anova.cca(dbrda_model, by = "axis", permutations = 999)
anova.cca(dbrda_model, by = "terms", permutations = 999)

# Extract scores and vectors
sample_scores <- as.data.frame(scores(dbrda_model, display = "sites"))
sample_scores$Sample <- rownames(sample_scores)
sample_scores <- left_join(sample_scores, df_combined %>% dplyr::select(Sample, Type, Treatment), by = "Sample")

phylum_vectors <- envfit(dbrda_model, df_phylum_hel, permutations = 999)
species_scores <- as.data.frame(phylum_vectors$vectors$arrows)
species_scores$phylum <- rownames(species_scores)
species_scores$length <- phylum_vectors$vectors$r
species_scores$CAP1 <- species_scores$CAP1 * species_scores$length
species_scores$CAP2 <- species_scores$CAP2 * species_scores$length

bp_scores <- as.data.frame(scores(dbrda_model, display = "bp"))
bp_scores$var <- rownames(bp_scores)

# Replace variable names for display
bp_scores$var <- recode(bp_scores$var,
                        "d13C" = "δ13C",
                        "d15N" = "δ15N")
# Define colors
type_colors <- c(
  "alder_l" = "#D73027",
  "alder_s" = "#FC8D59",
  "algae_l" = "#1A9850",
  "algae_s" = "#91CF60",
  "feces_l" = "#4575B4",
  "feces_s" = "#74ADD1"
)

# Visualization
p <- ggplot(sample_scores, aes(x = CAP1, y = CAP2, color = Type, shape = Treatment)) +
  geom_point(size = 5, alpha = 0.8) +
  geom_segment(data = bp_scores, 
               aes(x = 0, y = 0, xend = CAP1 * 1.2, yend = CAP2 * 1.2), 
               inherit.aes = FALSE, 
               arrow = arrow(length = unit(0.15, "cm")), 
               color = "black", 
               linewidth = 0.7) +
  geom_text_repel(data = bp_scores, 
                  aes(x = CAP1 * 1.4, y = CAP2 * 1.4, label = var),
                  inherit.aes = FALSE, 
                  color = "black", 
                  size = 5) +
  geom_segment(data = species_scores, 
               aes(x = 0, y = 0, xend = CAP1 * 1.2, yend = CAP2 * 1.2, linewidth = length), 
               inherit.aes = FALSE, 
               arrow = arrow(length = unit(0.1, "cm")), 
               color = "#F8766D", 
               linetype = "dashed") +
  geom_text_repel(data = species_scores, 
                  aes(x = CAP1 * 1.5, y = CAP2 * 1.5, label = phylum),
                  inherit.aes = FALSE, 
                  color = "#F8766D", 
                  size = 5,
                  box.padding = 0.5, 
                  point.padding = 0.5, 
                  segment.color = "grey50") +
  scale_linewidth_continuous(range = c(0.2, 0.6)) +

  scale_color_manual(values = type_colors) +
  theme_minimal(base_size = 16) +
  theme(
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 18, face = "bold")
  ) +
  labs(
    title = "dbRDA: Phylum vs Environmental Vectors (Bray-Curtis)(p < 0.05)",
    x = "CAP1", 
    y = "CAP2",
    color = "Type",
    shape = "Treatment",
    linewidth = "Impact (R²)"
  )

print(p)