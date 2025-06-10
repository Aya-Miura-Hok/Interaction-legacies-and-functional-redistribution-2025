# Load required libraries
library(phyloseq)
library(DESeq2)
library(ggplot2)
library(metagenomeSeq)  # CSS正規化用
library(ggrepel)  # ラベルの重なりを防ぐために使用
library(dplyr)

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
fpomsize <- subset_samples(exp, Treatment %in% c("C", "I", "T"))
exp_filtered <- prune_samples(sample_sums(fpomsize) > 1474, fpomsize)
otu_table_filtered <- prune_taxa(taxa_sums(exp_filtered) > (0.01 * nsamples(exp_filtered)), exp_filtered)
phylum_table <- tax_glom(otu_table_filtered, taxrank = "Phylum")

# Rename NA phylum
tax_data <- as.data.frame(tax_table(phylum_table))
otu_names <- taxa_names(phylum_table)
na_count <- 1
new_otu_names <- sapply(otu_names, function(otu) {
  phylum_name <- tax_data$Phylum[match(otu, rownames(tax_data))]
  if (is.na(phylum_name) || phylum_name == "") {
    new_name <- paste0("NA", na_count)
    na_count <<- na_count + 1
    return(new_name)
  } else {
    return(phylum_name)
  }
})
taxa_names(phylum_table) <- new_otu_names

# Comparison targets
comparison_list <- list(
  list(field = "Type2", group1 = "alder", group2 = "algae"),
  list(field = "Type2", group1 = "alder", group2 = "feces"),
  list(field = "Type2", group1 = "algae", group2 = "feces"),
  list(field = "Treatment", group1 = "C", group2 = "I"),
  list(field = "Treatment", group1 = "C", group2 = "T")
)

# Main loop for comparisons
for (cmp in comparison_list) {
  field <- cmp$field
  group1 <- cmp$group1
  group2 <- cmp$group2
  contrast_label <- paste(field, group1, "vs", group2)

  message("Running DESeq2 for: ", contrast_label)

  dds <- phyloseq_to_deseq2(phylum_table, as.formula(paste("~", field)))
  geoMeans <- apply(counts(dds), 1, function(x) exp(mean(log(x[x > 0]), na.rm=TRUE)))
  dds <- estimateSizeFactors(dds, geoMeans = geoMeans)
  dds <- DESeq(dds)

  res <- results(dds, contrast = c(field, group1, group2))
  res <- na.omit(res[order(res$padj), ])
  res_df <- as.data.frame(res)
  res_df$Phylum <- rownames(res_df)
  sig_res_df <- subset(res_df, padj < 0.05)

  p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = padj < 0.05, size = abs(log2FoldChange))) +
    geom_text_repel(data = sig_res_df, aes(label = Phylum), size = 4, box.padding = 0.5, max.overlaps = 15) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
    labs(
      title = paste0("Volcano Plot: ", group1, " vs ", group2, " (", field, ")"),
      x = "Log2 Fold Change",
      y = "-Log10 Adjusted p-value"
    ) +
    theme_minimal()

  print(p)
}