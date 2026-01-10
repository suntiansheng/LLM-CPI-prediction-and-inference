suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(lubridate)
  library(readr)
})

source("utility.R")
set.seed(1234)

df_raw <- read.csv("wsj_mean_inflation_predictions_withBERT_vector.csv",
                   check.names = FALSE, stringsAsFactors = FALSE)


target_var <- "CPIMEDSL"
macro_raw <- read.csv("./2025-10-MD.csv", stringsAsFactors = FALSE, check.names = FALSE)
macro_raw$sasdate <- suppressWarnings(as.Date(macro_raw[[1]], format = "%m/%d/%Y"))
macro_raw <- macro_raw[!is.na(macro_raw$sasdate), ]


topic_cols <- names(df_raw)[grepl("^cls_\\d+$", names(df_raw))]
if (length(topic_cols) == 0) {
  stop("No topic columns found in wsj_mean_inflation_predictions_withBERT_vector.csv")
}


df <- df_raw %>%
  mutate(date = as.Date(date)) %>%
  filter(date < as.Date("2025-09-01"), date > as.Date("2015-08-31")) %>%
  arrange(date) %>%
  mutate(
    month = format(date, "%Y-%m"),
    day = day(date),
    period = case_when(
      day <= 10 ~ "up",
      day <= 20 ~ "middle",
      TRUE ~ "down"
    )
  )

df_summary <- df %>%
  group_by(month, period) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

y_m_list <- sort(unique(df_summary$month))
month_dates <- as.Date(paste0(y_m_list, "-01"))


X_m <- matrix(ncol = length(topic_cols), nrow = length(y_m_list))
for (i in seq_along(y_m_list)) {
  res_ym <- y_m_list[i]
  X_m[i, ] <- apply(df_summary[df_summary$month == res_ym, topic_cols], 2, mean)
}
X_m <- scale(X_m)
colnames(X_m) <- topic_cols

topic_dates <- df$date
topic_months <- df$month

m_dates <- macro_raw$sasdate
if (!target_var %in% names(macro_raw)) {
  stop(paste("Target variable not found in 2025-10-MD.csv:", target_var))
}
m_vals <- macro_raw[[target_var]]
ptr <- 1
current <- NA_real_
target_series <- numeric(length(topic_dates))
for (i in seq_along(topic_dates)) {
  t_d <- topic_dates[i]
  while (ptr <= length(m_dates) && m_dates[ptr] <= t_d) {
    current <- m_vals[ptr]
    ptr <- ptr + 1
  }
  target_series[i] <- current
}

y_df <- data.frame(month = topic_months, y = target_series) %>%
  group_by(month) %>%
  summarise(y = mean(y, na.rm = TRUE), .groups = "drop") %>%
  filter(month %in% y_m_list)
y_raw <- y_df$y
y_ratio <- c(1, y_raw[-1] / y_raw[-length(y_raw)])
y_center <- mean(y_ratio, na.rm = TRUE)
y_scale <- sd(y_ratio, na.rm = TRUE)
y <- (y_ratio - y_center) / y_scale


EOD <- max(which(month_dates < as.Date("2024-06-01")))
train_idx <- 1:EOD


sel <- tryCatch(
  select_topics(y, X_m, train_idx = train_idx, use_cor = FALSE, penalty_w = 10),
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
  geom_jitter(width = 0.2, alpha = 0.4, size = 1.8, color = "grey60") +
  geom_point(data = subset(score_df, selected), aes(y = abs_cor), color = "#d95f02", size = 3) +
  labs(
    x = NULL,
    y = "|Correlation| with AR residuals",
    title = NULL
  ) +
  theme_minimal() +
  theme(
    legend.position = "top",
    plot.title = element_blank()
  )+coord_flip()

p

ggsave("boxplot.png", p, width = 7, height = 5, dpi = 300)
