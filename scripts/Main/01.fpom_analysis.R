# Load required libraries
library(phyloseq)
library(ggplot2)
library(dplyr)
library(car)
library(FSA)
library(rstatix)
library(rcompanion)
library(PMCMRplus)
library(multcomp)
library(reshape2)

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
ps_subset <- subset_samples(ps, Treatment %in% c("C", "I", "M"))

# Define FPOM types
fpom_types <- list(
  alder = c("alder_l", "alder_s"),
  algae = c("algae_l", "algae_s"),
  feces = c("feces_l", "feces_s")
)

# Main analysis function
analyze_fpom_auto <- function(ps_type, type_label) {
  df <- as.data.frame(as.matrix(sample_data(ps_type)))
  df$d13C <- as.numeric(as.character(df$d13C))
  df$CN <- as.numeric(as.character(df$CN))
  df$d15N <- as.numeric(as.character(df$d15N))
  df$Treatment <- factor(df$Treatment, levels = c("C", "I", "M"))
  df$Size <- ifelse(grepl("_l$", df$Type), "Large", "Small")
  df$Size <- as.factor(df$Size)

  cat("\n=== Analysis for:", toupper(type_label), "===\n")

  for (var in c("d13C", "CN", "d15N")) {
    cat("\n--- Variable:", var, "---\n")
    response <- df[[var]]

    # Normality and homogeneity
    normal_p <- shapiro.test(response)$p.value
    levene_p <- leveneTest(response ~ Treatment * Size, data = df)[1, "Pr(>F)"]

    cat("Shapiro-Wilk p =", round(normal_p, 4), "| Levene’s p =", round(levene_p, 4), "\n")

    is_normal <- normal_p > 0.05
    is_homoscedastic <- levene_p > 0.05

    if (is_normal && is_homoscedastic) {
      cat("→ Using parametric 2-way ANOVA + Dunnett's test\n")
      model <- aov(response ~ Treatment * Size, data = df)
      print(summary(model))

      posthoc <- glht(model, linfct = mcp(Treatment = "Dunnett"))
      print(summary(posthoc))

      for (sz in c("Large", "Small")) {
        df_sub <- df[df$Size == sz, ]
        cat(paste0("\nDunnett test for ", sz, " samples:\n"))
        model_sub <- aov(as.formula(paste(var, "~ Treatment")), data = df_sub)
        print(summary(glht(model_sub, linfct = mcp(Treatment = "Dunnett"))))
      }

    } else {
      cat("→ Using Scheirer–Ray–Hare + Dunn’s test (non-parametric Dunnett)\n")
      print(scheirerRayHare(as.formula(paste(var, "~ Treatment + Size + Treatment:Size")), data = df))

      for (sz in c("Large", "Small")) {
        df_sub <- df[df$Size == sz, ]
        cat(paste0("\nDunn’s test for ", sz, " samples:\n"))
        print(dunnettTest(x = df_sub[[var]], g = df_sub$Treatment, p.adjust.method = "dunnett"))
      }
    }
  }

  # Plotting
  df_long <- melt(df, id.vars = c("Treatment", "Size"),
                  measure.vars = c("d13C", "CN", "d15N"),
                  variable.name = "Variable", value.name = "Value")

  df_long$Treatment <- factor(df_long$Treatment, levels = c("C", "I", "M"))
  df_long$Variable <- factor(df_long$Variable, levels = c("CN", "d13C", "d15N"))
  levels(df_long$Variable) <- c("C:N", expression(delta^13*C), expression(delta^15*N))

  p <- ggplot(df_long, aes(x = Treatment, y = Value, fill = Size)) +
    geom_boxplot() +
    facet_wrap(~Variable, nrow = 1, scales = "free_y", labeller = label_parsed) +
    labs(title = paste("FPOM source:", type_label), x = "Treatment", y = "Value") +
    theme_bw() +
    theme(
      panel.border = element_rect(color = "black", fill = NA),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "gray80"),
      axis.text.x = element_text(size = 16),
      axis.text.y = element_text(size = 14),
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 11),
      strip.text = element_text(size = 13, face = "bold")
    )
  print(p)
}

# Run for each FPOM type
for (type_label in names(fpom_types)) {
  ps_type <- subset_samples(ps_subset, Type %in% fpom_types[[type_label]])
  analyze_fpom_auto(ps_type, type_label)
}
