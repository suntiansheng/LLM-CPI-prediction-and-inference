library(dplyr)
library(forecast)
library(MTS)
library(lubridate)
library(ggplot2)


data_path <- "topics_with_pred_score_final.csv"
df_raw <- read.csv(data_path, stringsAsFactors = FALSE, check.names = FALSE)

df <- df_raw %>%
  mutate(date = as.Date(date)) %>%
  arrange(date)


df_period <- df %>%
  mutate(
    month = format(date, "%Y-%m"),
    day = day(date),
    period = case_when(
      day <= 10 ~ "up",
      day <= 20 ~ "middle",
      TRUE ~ "down"
    )
  ) %>%
  dplyr::select(-date)

df_summary <- df_period %>%
  group_by(month, period) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

period_levels <- c("up", "middle", "down")
complete_months <- df_summary %>%
  count(month) %>%
  filter(n == length(period_levels)) %>%
  pull(month)
y_m_list <- sort(complete_months)
df_summary <- df_summary %>% filter(month %in% complete_months)


tilde_y_m <- matrix(ncol = length(period_levels), nrow = length(y_m_list))
for (i in seq_along(y_m_list)) {
  res_ym <- y_m_list[i]
  rows <- df_summary %>%
    filter(month == res_ym) %>%
    arrange(factor(period, levels = period_levels))
  if (nrow(rows) != length(period_levels)) {
    stop(sprintf("Month %s missing period rows; found %s", res_ym, nrow(rows)))
  }
  tilde_y_m[i, ] <- rows$pred_score
}


topic_cols <- names(df_summary)[grepl("^topic_\\d+$", names(df_summary))]
if (length(topic_cols) == 0) {
  stop("Topic columns missing from data.")
}

X_m <- matrix(ncol = length(topic_cols), nrow = length(y_m_list))
for (i in seq_along(y_m_list)) {
  res_ym <- y_m_list[i]
  X_m[i, ] <- df_summary %>%
    filter(month == res_ym) %>%
    dplyr::select(all_of(topic_cols)) %>%
    summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>%
    unlist(use.names = FALSE)
}


y_df <- df_period %>%
  group_by(month) %>%
  summarise(y = mean(CPIAUCSL, na.rm = TRUE), .groups = "drop") %>%
  filter(!is.na(y)) %>%
  arrange(month)


common_months <- intersect(y_m_list, y_df$month)
month_filter <- match(common_months, y_m_list)
tilde_y_m <- tilde_y_m[month_filter, , drop = FALSE]
X_m <- X_m[month_filter, , drop = FALSE]
y_df <- y_df %>% filter(month %in% common_months) %>%
  arrange(match(month, common_months))
y_m_list <- common_months

y <- y_df$y
y <- as.numeric(scale(y))

h <- 5
test_idx <- (length(y) - h + 1):length(y)
obs_idx <- setdiff(seq_along(y), test_idx)

ar_fit <- auto.arima(
  y[obs_idx],
  max.q = 0,
  D = 0,
  seasonal = FALSE,
  allowmean = FALSE,
  allowdrift = FALSE
)
p_use <- max(1, ar_fit$arma[1])
ar_fit <- arima(y[obs_idx], order = c(p_use, 0, 0), include.mean = FALSE)

strong_idx <- order(abs(cov(ar_fit$residuals, X_m[obs_idx, ])), decreasing = TRUE)[1:2]

arx_fit <- tryCatch(
  arima(y[obs_idx], order = c(p_use, 0, 0), xreg = X_m[obs_idx, strong_idx], include.mean = FALSE),
  error = function(e) {
    arima(y[obs_idx], order = c(1, 0, 0), xreg = X_m[obs_idx, strong_idx], include.mean = FALSE)
  }
)
arx_predictions <- predict(arx_fit, n.ahead = h, newxreg = X_m[test_idx, strong_idx])$pred

tilde_y_fit <- VARX(zt = tilde_y_m, p = 1, xt = X_m[, strong_idx], m = 0, include.mean = FALSE)
tilde_y_residual_m <- tilde_y_fit$residuals[obs_idx, ]

combined_data <- bind_rows(
  data.frame(x = as.numeric(arx_fit$residuals), y = as.numeric(tilde_y_residual_m[, 1]), group = "First"),
  data.frame(x = as.numeric(arx_fit$residuals), y = as.numeric(tilde_y_residual_m[, 2]), group = "Middle"),
  data.frame(x = as.numeric(arx_fit$residuals), y = as.numeric(tilde_y_residual_m[, 3]), group = "Last")
)

combined_data$group <- factor(combined_data$group, levels = c("First", "Middle", "Last"))

remove_outliers <- function(data, x_col, y_col, prob = 0.975) {
  coords <- data %>%
    dplyr::select(all_of(c(x_col, y_col))) %>%
    na.omit()
  if (nrow(coords) == 0) {
    return(data)
  }
  center <- colMeans(coords)
  cov_mat <- cov(coords)

  cov_mat <- cov_mat + diag(1e-6, nrow = 2)
  dist2 <- mahalanobis(coords, center = center, cov = cov_mat)
  cutoff <- qchisq(prob, df = 2)
  data %>%
    mutate(keep_flag = dist2 <= cutoff) %>%
    filter(keep_flag) %>%
    dplyr::select(-keep_flag)
}

combined_data_cleaned <- combined_data %>%
  group_by(group) %>%
  group_modify(~ remove_outliers(.x, "x", "y"))

band_x <- max(bw.nrd(combined_data_cleaned$x), 0.1)
band_y <- max(bw.nrd(combined_data_cleaned$y), 0.1)
band_width <- c(band_x * 3, band_y * 3)

p <- ggplot(combined_data_cleaned, aes(x = x, y = y)) +
  stat_density_2d(aes(color = after_stat(level)), linewidth = 1, h = band_width) +
  geom_point(alpha = 0.3, color = "blue") +
  scale_color_viridis_c() +
  facet_wrap(~group, ncol = 3) +
  labs(x = "LLM inflation index", y = "CPI", color = "Density level") +
  xlim(-.3, .3) + ylim(-.4, .4) +
  theme_minimal()
p

ggsave("error_density_plot.png", p, width = 10, height = 4, dpi = 300)
