# Load required libraries
library(phyloseq)
library(microbiome)
library(tidyverse)
library(ggpubr)

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

# Aggregate rare taxa and normalize to relative abundance (Phylum level)
pseq <- fpomsize %>%
  aggregate_rare(level = "Phylum", detection = 5, prevalence = 0.4) %>%
  microbiome::transform("compositional") %>%
  aggregate_taxa(level = "Phylum")

# Subset Verrucomicrobiota
Verru <- subset_taxa(pseq, Phylum == "Verrucomicrobiota")

# Extract environmental variables
df_env <- as(sample_data(Verru), "data.frame") %>%
  rownames_to_column(var = "Sample") %>%
  select(Sample, d13C, HIX, BIX)

# Extract relative abundance of Verrucomicrobiota
verru_df <- psmelt(Verru) %>%
  group_by(Sample) %>%
  summarise(Verrucomicrobiota = sum(Abundance)) %>%
  ungroup()

# Merge with environmental variables
df_lm <- df_env %>%
  left_join(verru_df, by = "Sample") %>%
  drop_na()

# Linear models
model_d13C <- lm(d13C ~ Verrucomicrobiota, data = df_lm)
model_HIX  <- lm(HIX  ~ Verrucomicrobiota, data = df_lm)
model_BIX  <- lm(BIX  ~ Verrucomicrobiota, data = df_lm)

# Model summaries
summary(model_d13C)
summary(model_HIX)
summary(model_BIX)

# Convert to long format for visualization
df_long <- df_lm %>%
  pivot_longer(cols = c(d13C, HIX, BIX), names_to = "Variable", values_to = "Value") %>%
  mutate(Variable = factor(Variable, levels = c("d13C", "HIX", "BIX")))

# Plot: scatter plots with regression lines and correlation statistics
ggplot(df_long, aes(x = Verrucomicrobiota, y = Value)) +
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
    title = "Linear relationships between Verrucomicrobiota and DOM metrics"
  )