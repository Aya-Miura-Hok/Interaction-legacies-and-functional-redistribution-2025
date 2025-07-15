# Load required packages
library(dplyr)
library(tidyr)

# Load DOC concentration data
df <- read.csv("https://raw.githubusercontent.com/Aya-Miura-Hok/Interaction-legacies-and-functional-redistribution-2025/refs/heads/main/data/DOC_conc.csv", stringsAsFactors = FALSE)

# Standardize DOC column name
colnames(df)[which(grepl("DOC", colnames(df)))] <- "DOC"

# Calculate group-wise mean and standard deviation
summary_stats <- df %>%
  group_by(Type, Size) %>%
  summarise(
    mean = round(mean(DOC), 1),
    sd = round(sd(DOC), 1),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = Size, values_from = c(mean, sd))

# Perform t-tests for each Type (Large vs Small)
t_test_results <- df %>%
  group_split(Type) %>%
  lapply(function(group_data) {
    test <- t.test(DOC ~ Size, data = group_data)
    data.frame(
      Type = unique(group_data$Type),
      p_value = if (test$p.value < 0.01) {
        "p < 0.01"
      } else {
        paste0("p = ", formatC(test$p.value, digits = 2, format = "f"))
      }
    )
  }) %>%
  bind_rows()

# Merge summary statistics with t-test results
summary_table <- left_join(summary_stats, t_test_results, by = "Type") %>%
  select(
    Type,
    DOC_L = mean_L, DOC_S = mean_S,
    SD_L = sd_L, SD_S = sd_S,
    p_value
  )

# Display table
print(summary_table)