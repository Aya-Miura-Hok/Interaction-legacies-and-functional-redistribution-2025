library(dplyr)
library(tidyr)

# データ読み込み
df <- read.csv("https://raw.githubusercontent.com/Aya-Miura-Hok/Interaction-legacies-and-functional-redistribution-2025/refs/heads/main/data/DOC_conc.csv", stringsAsFactors = FALSE)

# DOC列名を明示
colnames(df)[which(grepl("DOC", colnames(df)))] <- "DOC"

# 平均・SDを計算
summary_stats <- df %>%
  group_by(Type, Size) %>%
  summarise(mean = round(mean(DOC), 1),
            sd = round(sd(DOC), 1),
            .groups = "drop") %>%
  pivot_wider(names_from = Size, values_from = c(mean, sd))

# Typeごとのt検定結果を格納
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

# 平均・SDとp値を結合
summary_table <- left_join(summary_stats, t_test_results, by = "Type") %>%
  select(Type,
         DOC_L = mean_L, DOC_S = mean_S,
         SD_L = sd_L, SD_S = sd_S,
         p_value)

# 表示
print(summary_table)