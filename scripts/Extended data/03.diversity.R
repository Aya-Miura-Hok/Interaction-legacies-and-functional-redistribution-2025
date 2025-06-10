# Load required libraries
library(phyloseq)
library(iNEXT)
library(dplyr)
library(ggplot2)
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
exp<-subset_samples(ps, Type %in% c("alder_l", "alder_s", "algae_l", "algae_s", "feces_l", "feces_s"))
ps_subset <- subset_samples(exp, Treatment %in% c("C", "I", "T"))
physeq_filtered <- prune_samples(sample_sums(ps_subset) > 999, ps_subset)

# Map Type names to integers
type_names <- unique(sample_data(physeq_filtered)$Type)
type_mapping <- setNames(seq_along(type_names), type_names)
sample_data(physeq_filtered)$Type <- unlist(lapply(sample_data(physeq_filtered)$Type, function(x) type_mapping[x]))

# iNEXT analysis per Type
type_list <- list()
inext_results <- list()

for (i in seq_along(type_mapping)) {
  physeq_subset <- subset_samples(physeq_filtered, Type == i)
  otu_mat <- as.data.frame(otu_table(physeq_subset))
  otu_trans <- t(otu_mat)
  sample_counts <- apply(otu_trans, 2, function(x) x[x > 0])
  inext_results[[as.character(i)]] <- iNEXT(sample_counts, q = c(0, 1, 2), datatype = "abundance")
  type_list[[as.character(i)]] <- physeq_subset
}

# Plot diversity curves
par(mfrow = c(2, 2))
for (i in seq_along(inext_results)) {
  plot(inext_results[[i]], type = 1, main = paste("Type:", names(type_mapping)[i]))
}

# Aggregate iNEXT results
type_labels <- names(type_mapping)
all_results_size_based <- list()
all_results_coverage_based <- list()
all_colnames_size <- character()
all_colnames_coverage <- character()

for (i in seq_along(type_labels)) {
  res <- inext_results[[as.character(i)]]
  type <- type_labels[i]
  if (!is.null(res$iNextEst)) {
    if (!is.null(res$iNextEst$size_based)) {
      df_size <- as.data.frame(res$iNextEst$size_based)
      df_size$Type <- type
      all_colnames_size <- union(all_colnames_size, colnames(df_size))
      all_results_size_based[[length(all_results_size_based) + 1]] <- df_size
    }
    if (!is.null(res$iNextEst$coverage_based)) {
      df_cov <- as.data.frame(res$iNextEst$coverage_based)
      df_cov$Type <- type
      all_colnames_coverage <- union(all_colnames_coverage, colnames(df_cov))
      all_results_coverage_based[[length(all_results_coverage_based) + 1]] <- df_cov
    }
  }
}

# Standardize column structure and bind
for (i in seq_along(all_results_size_based)) {
  missing <- setdiff(all_colnames_size, names(all_results_size_based[[i]]))
  for (col in missing) all_results_size_based[[i]][[col]] <- NA
  all_results_size_based[[i]] <- all_results_size_based[[i]][, all_colnames_size, drop = FALSE]
}

for (i in seq_along(all_results_coverage_based)) {
  missing <- setdiff(all_colnames_coverage, names(all_results_coverage_based[[i]]))
  for (col in missing) all_results_coverage_based[[i]][[col]] <- NA
  all_results_coverage_based[[i]] <- all_results_coverage_based[[i]][, all_colnames_coverage, drop = FALSE]
}

# Output combined coverage-based result
final_coverage_df <- bind_rows(all_results_coverage_based)
write.csv(final_coverage_df, "iNEXT_Sitewise_Results_coverage_based_type.csv", row.names = FALSE)

# Load required libraries
library(dplyr)
library(ggplot2)
library(FSA)
library(ggsignif)
library(tidyr)
library(purrr)
library(gridExtra)

# Load coverage-based diversity estimates
coverage_based_results <- read.csv("https://raw.githubusercontent.com/Aya-Miura-Hok/Aya-Miura-Hok/Interaction-legacies-and-functional-redistribution-2025/refs/heads/main/data/bacteria/iNEXT_Sitewise_Results_coverage_based_type.csv")

# Filter rows close to SC = 0.95
coverage_95_results <- dplyr::filter(coverage_based_results, abs(SC - 0.95) < 0.01)

# Extract diversity estimates for q = 0, 1, 2
diversity_at_95_coverage <- coverage_95_results %>%
  dplyr::filter(Order.q %in% c(0, 1, 2)) %>%
  dplyr::select(Type, Order.q, qD, qD.LCL, qD.UCL)

# Initialize plot container
plots <- list()

# Loop through each q-order
for (q in c(0, 1, 2)) {
  data_q <- dplyr::filter(diversity_at_95_coverage, Order.q == q)
  data_q$Type <- factor(data_q$Type)

  kruskal_test <- kruskal.test(qD ~ Type, data = data_q)
  dunn_test <- dunnTest(qD ~ Type, data = data_q, method = "bonferroni")

  # Handle varying column name conventions
  if ("Comparison" %in% colnames(dunn_test$res) & "P.adj" %in% colnames(dunn_test$res)) {
    significant_pairs <- dunn_test$res %>%
      dplyr::filter(P.adj < 0.05) %>%
      dplyr::select(Comparison, P.adj)
    p_col <- "P.adj"
    name_col <- "Comparison"
  } else {
    significant_pairs <- dunn_test$res %>%
      dplyr::filter(p.adj < 0.05) %>%
      dplyr::select(comparison, p.adj)
    p_col <- "p.adj"
    name_col <- "comparison"
  }

  if (nrow(significant_pairs) > 0) {
    signif_annotations <- significant_pairs %>%
      tidyr::separate(1, into = c("Type1", "Type2"), sep = " - ") %>%
      mutate(
        y_position = max(data_q$qD, na.rm = TRUE) + (1:n()) * 0.1 * max(data_q$qD, na.rm = TRUE),
        annotations = if_else(.[[p_col]] < 0.01, "**", "*")
      )

    comparisons_list <- signif_annotations %>%
      dplyr::select(Type1, Type2) %>%
      split(1:nrow(.)) %>%
      purrr::map(~ as.character(c(.x$Type1, .x$Type2)))

    plot <- ggplot(data_q, aes(x = Type, y = qD, fill = Type)) +
      geom_boxplot() +
      labs(
        title = paste0("Order.q = ", q),
        x = "Type",
        y = "qD (Diversity Index)"
      ) +
      theme_minimal(base_size = 16) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        legend.position = "none"
      ) +
      geom_signif(
        comparisons = comparisons_list,
        annotations = signif_annotations$annotations,
        y_position = signif_annotations$y_position,
        tip_length = 0.01,
        textsize = 6
      )
  } else {
    plot <- ggplot(data_q, aes(x = Type, y = qD, fill = Type)) +
      geom_boxplot() +
      labs(
        title = paste0("Order.q = ", q),
        x = "Type",
        y = "qD (Diversity Index)"
      ) +
      theme_minimal(base_size = 16) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        legend.position = "none"
      )
  }

  plots[[paste0("Order.q_", q)]] <- plot
}

# Display plots side by side
p <- gridExtra::grid.arrange(
  plots$Order.q_0,
  plots$Order.q_1,
  plots$Order.q_2,
  nrow = 1
)

print(p)