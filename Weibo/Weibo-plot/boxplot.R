suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(lubridate)
  library(readr)
})

source("utility.R")
set.seed(1234)

df_raw <- read.csv("bert_final.csv", check.names = FALSE, stringsAsFactors = FALSE)


topic_cols <- names(df_raw)[grepl("^X\\d+$|^\\d+$", names(df_raw))]
if (length(topic_cols) == 0) {
  stop("No topic columns found in bert_final.csv")
}


df_new <- df_raw %>%
  mutate(
    time = as.Date(time),
    month = format(time, "%Y-%m"),
    day = day(time),
    period = case_when(
      day <= 10 ~ "up",
      day <= 20 ~ "middle",
      TRUE ~ "down"
    ),
    unem = unem,
    cpi = cpi_lastmonth_100
  ) %>%
  dplyr::select(month, period, aveScore, all_of(topic_cols), unem, cpi)

df_summary <- df_new %>%
  group_by(month, period) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

y_m_list <- sort(unique(df_summary$month))
month_dates <- as.Date(paste0(y_m_list, "-01"))


topic_cols_summary <- intersect(topic_cols, names(df_summary))
X_m <- df_summary %>%
  group_by(month) %>%
  summarise(across(all_of(topic_cols_summary), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  arrange(month) %>%
  dplyr::select(-month) %>%
  as.matrix()
X_m <- scale(X_m)


y <- df_new %>%
  group_by(month) %>%
  summarise(y = mean(cpi, na.rm = TRUE), .groups = "drop") %>%
  arrange(month) %>%
  pull(y)
y <- scale(y - 100)


EOD <- max(which(month_dates < as.Date("2024-07-01")))
train_idx <- 1:EOD


sel <- tryCatch(
  select_topics(y, X_m, train_idx = train_idx, use_cor = FALSE, penalty_w = 0.1),
  error = function(e) list(selected = character())
)
selected_topics <- sel$selected


ar_base <- auto.arima(y[train_idx], max.q = 0, D = 0, seasonal = FALSE,
                      allowmean = FALSE, allowdrift = FALSE)
p_ar <- ar_base$arma[1]
ar_fit <- arima(y[train_idx], order = c(p_ar, 0, 0), include.mean = FALSE)
resid_base <- residuals(ar_fit)

cor_scores <- sapply(colnames(X_m), function(col) {
  cor(resid_base, X_m[train_idx, col], use = "complete.obs")
})

score_df <- data.frame(
  topic = colnames(X_m),
  abs_cor = abs(cor_scores),
  selected = colnames(X_m) %in% selected_topics,
  stringsAsFactors = FALSE
)


p <- ggplot(score_df, aes(x = "All topics", y = abs_cor)) +
  geom_boxplot(outlier.shape = NA, fill = "grey90", color = "grey40") +
  geom_jitter(aes(color = selected), width = 0.2, alpha = 0.4, size = 1.8) +
  geom_point(data = subset(score_df, selected), aes(y = abs_cor), color = "#d95f02", size = 3) +
  labs(
    x = NULL,
    y = "|Correlation| with AR residuals",
    color = "LLM-TS selected",
    title = NULL
  ) +
  scale_color_manual(values = c(`TRUE` = "#d95f02", `FALSE` = "grey50")) +
  theme_minimal() +
  theme(
    legend.position = "top",
    plot.title = element_blank()
  )+coord_flip()

p

ggsave("boxplot.png", p, width = 7, height = 5, dpi = 300)
