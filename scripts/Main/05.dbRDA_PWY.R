# Load necessary libraries
library(phyloseq)
library(ggplot2)
library(dplyr)
library(microbiome)
library(vegan)
library(ggrepel)
library(tibble)
library(tidyr)

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

# Select key phyla from GLM results
selected_phyla <- c("Proteobacteria", "Verrucomicrobiota", "Planctomycetota", "Desulfobacterota", "Bacteroidota", "Cyanobacteria")
df_phylum_filtered <- dplyr::select(otu_table_long, Sample, all_of(selected_phyla))

# Load KO and pathway abundance tables
df_ko <- read.csv("https://raw.githubusercontent.com/Aya-Miura-Hok/Interaction-legacies-and-functional-redistribution-2025/refs/heads/main/data/bacteria/ko-feature-table.biom.csv", row.names = 1, check.names = FALSE)
df_pathabun <- read.csv("https://raw.githubusercontent.com/Aya-Miura-Hok/Interaction-legacies-and-functional-redistribution-2025/refs/heads/main/data/bacteria/pathabunｰfeature-table.biom.csv", row.names = 1, check.names = FALSE)

# Normalize each sample to sum to 1
df_ko_norm <- df_ko / rowSums(df_ko)
df_pathabun_norm <- df_pathabun / rowSums(df_pathabun)
df_ko_norm$Sample <- rownames(df_ko_norm)
df_pathabun_norm$Sample <- rownames(df_pathabun_norm)

# Select relevant KO and pathway IDs
selected_kos <- c("K01555", "K00480", "K01183", "K01225", "K02559", "K01785", "K00022", "K01067", "K01799", "K01301", "K01476", "K01206")
selected_pathways <- c("PWY-6590", "PWY-6628", "PWY-5453", "PWY-7398", "PWY-7219", "PWY-6386", "PWY-7074", "PWY-5386", "PWY-7299")

valid_kos <- intersect(selected_kos, colnames(df_ko_norm))
valid_pathways <- intersect(selected_pathways, colnames(df_pathabun_norm))

# Subset and join KO and pathway data
df_ko_filtered <- df_ko_norm %>% dplyr::select(Sample, all_of(valid_kos))
df_pathabun_filtered <- df_pathabun_norm %>% dplyr::select(Sample, all_of(valid_pathways))
df_functional <- inner_join(df_pathabun_filtered, df_ko_filtered, by = "Sample")

# Combine all data
df_combined_all <- inner_join(df_phylum_filtered, df_functional, by = "Sample")
df_phylum <- df_combined_all %>% dplyr::select(Sample, all_of(selected_phyla)) %>% column_to_rownames(var = "Sample")
df_ko_pwy <- df_combined_all %>% dplyr::select(Sample, starts_with("K"), starts_with("PWY")) %>% column_to_rownames(var = "Sample")

# Hellinger transformation
df_ko_pwy[is.na(df_ko_pwy)] <- 0
df_ko_pwy[is.nan(as.matrix(df_ko_pwy))] <- 0
df_ko_pwy <- df_ko_pwy[rowSums(df_ko_pwy) > 0, ]
df_ko_pwy_hel <- decostand(df_ko_pwy, method = "hellinger")

# Match rows
df_phylum_matched <- df_phylum[rownames(df_ko_pwy_hel), ]

# Perform dbRDA
dbrda_model <- capscale(vegdist(df_ko_pwy_hel, method = "bray") ~ ., data = df_phylum_matched)
summary(dbrda_model)
anova.cca(dbrda_model, permutations = 999)
anova.cca(dbrda_model, by = "axis", permutations = 999)
anova.cca(dbrda_model, by = "terms", permutations = 999)

# Extract scores
species_scores <- as.data.frame(scores(envfit(dbrda_model, df_ko_pwy_hel, permutations = 999)$vectors$arrows))
species_scores$KO_PWY <- rownames(species_scores)
species_scores$length <- envfit(dbrda_model, df_ko_pwy_hel)$vectors$r
species_scores$CAP1 <- species_scores$CAP1 * species_scores$length
species_scores$CAP2 <- species_scores$CAP2 * species_scores$length

bp_scores <- as.data.frame(scores(dbrda_model, display = "bp"))
sample_scores <- data.frame(scores(dbrda_model, display = "sites"))
sample_scores$Sample <- rownames(sample_scores)

# Join sample info
df_env <- as(sample_data(fpomsize), "data.frame") %>%
  rownames_to_column("Sample") %>%
  dplyr::select(Sample, Type, Treatment, d13C, HIX, BIX)

df_combined_all <- left_join(df_combined_all, df_env, by = "Sample")

sample_scores <- left_join(
  sample_scores,
  df_combined_all %>% dplyr::select(Sample, Type, Treatment),
  by = "Sample"
)

# Define color and shape
type_colors <- c("alder_l" = "#D73027", "alder_s" = "#FC8D59",
                 "algae_l" = "#1A9850", "algae_s" = "#91CF60",
                 "feces_l" = "#4575B4", "feces_s" = "#74ADD1")

treatment_shapes <- c("C" = 21, "I" = 22, "M" = 24)

# Plot
p <- ggplot() +
  geom_point(data = sample_scores, aes(x = CAP1, y = CAP2, fill = Type, shape = Treatment),
             size = 3, alpha = 0.8, color = "black") +
  geom_segment(data = species_scores, aes(x = 0, y = 0, xend = CAP1 * 1.5, yend = CAP2 * 1.5),
               arrow = arrow(length = unit(0.1, "cm")), color = "#F8766D", linetype = "dashed") +
  geom_text_repel(data = species_scores, aes(x = CAP1 * 1.7, y = CAP2 * 1.7, label = KO_PWY),
                  color = "#F8766D", size = 5) +
  geom_segment(data = bp_scores, aes(x = 0, y = 0, xend = CAP1 * 1.2, yend = CAP2 * 1.2),
               arrow = arrow(length = unit(0.2, "cm")), color = "black", linetype = "solid") +
  geom_text_repel(data = bp_scores, aes(x = CAP1 * 1.3, y = CAP2 * 1.3, label = rownames(bp_scores)),
                  color = "black", size = 5) +
  scale_fill_manual(values = type_colors) +
  scale_shape_manual(values = treatment_shapes) +
  guides(
    fill = guide_legend(override.aes = list(shape = 21)),  
    shape = guide_legend()                                 
  ) +
  theme_minimal(base_size = 16) +
  theme(
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
  ) +
  labs(
    title = "dbRDA: Phylum vs KO/PWY (p < 0.05)",
    x = "CAP1",
    y = "CAP2",
    fill = "Type",
    shape = "Treatment"
  )
  
print(p)

# Extract significant phyla (p < 0.05)
anova_terms <- anova.cca(dbrda_model, by = "terms", permutations = 999)
sig_phyla <- rownames(anova_terms)[anova_terms$`Pr(>F)` < 0.05]
sig_phyla <- intersect(sig_phyla, colnames(df_phylum_matched))

# Subset data and re-run dbRDA model
if (length(sig_phyla) >= 1) {
  df_phylum_sig <- df_phylum_matched[, sig_phyla, drop = FALSE]
  dbrda_sig <- capscale(vegdist(df_ko_pwy_hel, method = "bray") ~ ., data = df_phylum_sig)

# Summary and permutation tests
cat("=== dbRDA Model (Significant Phyla Only) ===\n")
  print(summary(dbrda_sig))
  
  cat("\n=== Permutation Test (Global) ===\n")
  print(anova.cca(dbrda_sig, permutations = 999))
  
  cat("\n=== Permutation Test (by Axis) ===\n")
  print(anova.cca(dbrda_sig, by = "axis", permutations = 999))
  
  cat("\n=== Permutation Test (by Term) ===\n")
  print(anova.cca(dbrda_sig, by = "terms", permutations = 999))
  
} else {
  cat("No significant phyla detected (p < 0.05). Skipping filtered dbRDA model.\n")
}

# Subset phylum data to only significant ones
df_phylum_sig <- df_phylum_matched[, sig_phyla, drop = FALSE]

# Recalculate dbRDA with only significant phyla
dbrda_sig <- capscale(vegdist(df_ko_pwy_hel, method = "bray") ~ ., data = df_phylum_sig)

# Updated species and sample scores
species_scores <- as.data.frame(scores(envfit(dbrda_sig, df_ko_pwy_hel, permutations = 999)$vectors$arrows))
species_scores$KO_PWY <- rownames(species_scores)
species_scores$length <- envfit(dbrda_sig, df_ko_pwy_hel)$vectors$r
species_scores$CAP1 <- species_scores$CAP1 * species_scores$length
species_scores$CAP2 <- species_scores$CAP2 * species_scores$length

bp_scores <- as.data.frame(scores(dbrda_sig, display = "bp"))
sample_scores <- data.frame(scores(dbrda_sig, display = "sites"))
sample_scores$Sample <- rownames(sample_scores)

# Re-join sample metadata
sample_scores <- left_join(sample_scores,
                           df_combined_all %>% dplyr::select(Sample, Type, Treatment),
                           by = "Sample")

# Define color and shape
type_colors <- c("alder_l" = "#D73027", "alder_s" = "#FC8D59",
                 "algae_l" = "#1A9850", "algae_s" = "#91CF60",
                 "feces_l" = "#4575B4", "feces_s" = "#74ADD1")

treatment_shapes <- c("C" = 21, "I" = 22, "M" = 24)

# Plot
p <- ggplot() +
  geom_point(data = sample_scores, aes(x = CAP1, y = CAP2, fill = Type, shape = Treatment),
             size = 3, alpha = 0.8, color = "black") +
  geom_segment(data = species_scores, aes(x = 0, y = 0, xend = CAP1 * 1.5, yend = CAP2 * 1.5),
               arrow = arrow(length = unit(0.1, "cm")), color = "#F8766D", linetype = "dashed") +
  geom_text_repel(data = species_scores, aes(x = CAP1 * 1.7, y = CAP2 * 1.7, label = KO_PWY),
                  color = "#F8766D", size = 5) +
  geom_segment(data = bp_scores, aes(x = 0, y = 0, xend = CAP1 * 1.2, yend = CAP2 * 1.2),
               arrow = arrow(length = unit(0.2, "cm")), color = "black", linetype = "solid") +
  geom_text_repel(data = bp_scores, aes(x = CAP1 * 1.3, y = CAP2 * 1.3, label = rownames(bp_scores)),
                  color = "black", size = 5) +
  scale_fill_manual(values = type_colors) +
  scale_shape_manual(values = treatment_shapes) +
  guides(
    fill = guide_legend(override.aes = list(shape = 21)),
    shape = guide_legend()                               
  ) +
  theme_minimal(base_size = 16) +
  theme(
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
  ) +
  labs(
    title = "dbRDA: Phylum vs KO/PWY (p < 0.05)",
    x = "CAP1",
    y = "CAP2",
    fill = "Type",
    shape = "Treatment"
  )
  
print(p)