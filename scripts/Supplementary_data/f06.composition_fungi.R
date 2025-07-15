# Load required libraries
library(ggplot2)
library(microViz)
library(microbiome)
library(phyloseq)
library(RColorBrewer)

# Set default ggplot theme
theme_set(theme_bw())

# Import data
sample_sheet <- read.csv("https://raw.githubusercontent.com/Aya-Miura-Hok/Interaction-legacies-and-functional-redistribution-2025/refs/heads/main/data/fungi/metadata_table.csv", row.names = 1)
asv_sheet <- read.csv("https://raw.githubusercontent.com/Aya-Miura-Hok/Interaction-legacies-and-functional-redistribution-2025/refs/heads/main/data/fungi/otu_table.csv", row.names = 1)
tax_sheet <- read.csv("https://raw.githubusercontent.com/Aya-Miura-Hok/Interaction-legacies-and-functional-redistribution-2025/refs/heads/main/data/fungi/taxonomy_table.csv", row.names = 1)

# Create phyloseq object
ps_all <- phyloseq(
  otu_table(asv_sheet, taxa_are_rows = FALSE),
  sample_data(sample_sheet),
  tax_table(as.matrix(tax_sheet))
)

# Filter bacterial data and subset treatment samples
ps<-subset_taxa(ps_all, Kingdom == "Fungi")
exp<-subset_samples(ps, Type %in% c("alder_l", "alder_s", "algae_l", "algae_s", "feces_l", "feces_s"))
fpomsize <- subset_samples(exp, Treatment %in% c("C", "I", "M"))

# ---------------- Phylum-level community composition ---------------- #
# Aggregate rare taxa, normalize by compositional transformation, and merge by 'Description'
pseq_phylum <- fpomsize %>%
  aggregate_rare(level = "Phylum", detection = 5, prevalence = 0.4) %>%
  microbiome::transform("compositional") %>%
  aggregate_taxa("Phylum") %>%
  merge_samples("Description") %>%
  microbiome::transform("compositional")

# Define custom color palette
phylum_colors <- c(
  "#CAE0AB", "#7BAFDE", "#4EB265", "#F6C141",
  "#AA6F9E", "#DC050C", "#437DBF", "#CAACCB",
  "#F7F056", "#90C987", "#1965B0", "#D9CCE3",
  "#72190E", "#F1932D"
)

# Generate and display the plot
plot_phylum <- plot_composition(pseq_phylum, otu.sort = "abundance", transform = "identity") +
  scale_fill_manual(values = phylum_colors) +
  guides(fill = guide_legend(title = "Phylum")) +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 12)
  )

print(plot_phylum)

# ---------------- Genus-level community composition ---------------- #
# Aggregate rare genera, normalize, and merge by 'Description'
pseq_genus <- fpomsize %>%
  aggregate_rare(level = "Genus", detection = 5, prevalence = 0.2) %>%
  microbiome::transform("compositional") %>%
  aggregate_taxa("Genus") %>%
  merge_samples("Description") %>%
  microbiome::transform("compositional")

# Reuse the color palette (expand if necessary)
genus_colors <- phylum_colors

# Generate and display the plot
plot_genus <- plot_composition(pseq_genus, otu.sort = "abundance", transform = "identity") +
  scale_fill_manual(values = genus_colors) +
  guides(fill = guide_legend(title = "Genus")) +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 12)
  )

print(plot_genus)