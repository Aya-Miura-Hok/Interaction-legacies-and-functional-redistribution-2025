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

# Import data from GitHub
sample_sheet <- read.csv("https://raw.githubusercontent.com/Aya-Miura-Hok/Interaction-legacies-and-functional-redistribution-2025/refs/heads/main/data/bacteria/sample_sheet.csv", row.names = 1)
asv_sheet <- read.csv("https://raw.githubusercontent.com/Aya-Miura-Hok/Interaction-legacies-and-functional-redistribution-2025/refs/heads/main/data/bacteria/otu_table.csv", row.names = 1)
tax_sheet <- read.csv("https://raw.githubusercontent.com/Aya-Miura-Hok/Interaction-legacies-and-functional-redistribution-2025/refs/heads/main/data/bacteria/taxonomy_table.csv", row.names = 1)

# Build phyloseq object
ps_all <- phyloseq(
  otu_table(asv_sheet, taxa_are_rows = FALSE),
  sample_data(sample_sheet),
  tax_table(as.matrix(tax_sheet))
)

ps <- subset_taxa(ps_all, Kingdom == "d__Bacteria")
ps_subset <- subset_samples(ps, Treatment %in% c("C", "I", "M"))

fpom_types <- list(
  alder = c("alder_l", "alder_s"),
  algae = c("algae_l", "algae_s"),
  feces = c("feces_l", "feces_s")
)

# Define function
analyze_HIX_BIX <- function(ps_type, type_label) {
  df <- as.data.frame(as.matrix(sample_data(ps_type)))
  df$HIX <- as.numeric(as.character(df$HIX))
  df$BIX <- as.numeric(as.character(df$BIX))
  df$Treatment <- factor(df$Treatment, levels = c("C", "I", "M"))
  df$Size <- ifelse(grepl("_l$", df$Type), "Large", "Small")
  df$Size <- as.factor(df$Size)

  cat("\n=== ANALYSIS FOR:", toupper(type_label), "===\n")

  for (var in c("HIX", "BIX")) {
    cat("\n--- Variable:", var, "---\n")
    response <- df[[var]]

    # Normality & homoscedasticity tests
    cat("Shapiro-Wilk normality:\n"); print(shapiro.test(response))
    cat("Levene’s test for homogeneity:\n")
    print(leveneTest(response ~ Treatment * Size, data = df))

    # Choose parametric or non-parametric based on Shapiro & Levene
    is_normal <- shapiro.test(response)$p.value > 0.05
    is_homoscedastic <- leveneTest(response ~ Treatment * Size, data = df)[1,"Pr(>F)"] > 0.05

    if (is_normal & is_homoscedastic) {
      cat("→ Using 2-way ANOVA + Dunnett's test\n")

      model <- aov(response ~ Treatment * Size, data = df)
      print(summary(model))
      posthoc <- glht(model, linfct = mcp(Treatment = "Dunnett"))
      print(summary(posthoc))

      for (sz in c("Large", "Small")) {
        df_sub <- subset(df, Size == sz)
        sub_model <- aov(as.formula(paste(var, "~ Treatment")), data = df_sub)
        cat("\nDunnett Test for", sz, "only:\n")
        print(summary(glht(sub_model, linfct = mcp(Treatment = "Dunnett"))))
      }

    } else {
      cat("→ Using Scheirer–Ray–Hare + non-parametric Dunnett (Dunn’s test)\n")
      print(scheirerRayHare(as.formula(paste(var, "~ Treatment + Size + Treatment:Size")), data = df))

      for (sz in c("Large", "Small")) {
        df_sub <- subset(df, Size == sz)
        cat("\nDunn’s Test for", sz, "only:\n")
        print(dunnettTest(x = df_sub[[var]], g = df_sub$Treatment, p.adjust.method = "dunnett"))
      }
    }
  }

  # Plotting
  df_long <- melt(df, id.vars = c("Treatment", "Size"),
                  measure.vars = c("HIX", "BIX"),
                  variable.name = "Variable",
                  value.name = "Value")

  df_long$Treatment <- factor(df_long$Treatment, levels = c("C", "I", "M"))

  p <- ggplot(df_long, aes(x = Treatment, y = Value, fill = Size)) +
    geom_boxplot() +
    facet_wrap(~Variable, nrow = 1, scales = "free_y") +
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

# Run analysis for all FPOM types
for (type_label in names(fpom_types)) {
  ps_type <- subset_samples(ps_subset, Type %in% fpom_types[[type_label]])
  analyze_HIX_BIX(ps_type, type_label)
}