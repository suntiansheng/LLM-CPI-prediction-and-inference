
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(forecast)
  library(MTS)
  library(ggplot2)
  library(patchwork)
})

source("./utility.R")


df <- read.csv("./topics_with_pred_score_final.csv")
df <- df[(df$date < "2025-09-01") & ("2015-08-31" < df$date), ]
df$date <- as.Date(df$date)
df <- df %>%
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
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE)))

y_m_list <- sort(unique(df_summary$month))
tilde_y_m <- matrix(ncol = 3, nrow = length(y_m_list))
for (i in seq_along(y_m_list)) {
  res_ym <- y_m_list[i]
  tilde_y_m[i, ] <- df_summary[df_summary$month == res_ym, ]$pred_score
}

topic_cols <- names(df_summary)[grepl("^topic", names(df_summary))]
X_m <- matrix(ncol = length(topic_cols), nrow = length(y_m_list))
for (i in seq_along(y_m_list)) {
  res_ym <- y_m_list[i]
  X_m[i, ] <- apply(df_summary[df_summary$month == res_ym, topic_cols], 2, mean)
}
X_m <- scale(X_m)
colnames(X_m) <- topic_cols


make_ratio <- function(x) {
  if (length(x) == 0) return(numeric())
  c(1, x[-1] / x[-length(x)])
}

targets <- list(
  list(col = "WPSFD49207", label = "PPI: Finished Goods"),
  list(col = "WPSFD49502", label = "PPI: Intermediate Materials"),
  list(col = "PPICMM", label = "PPI: Commodities"),
  list(col = "CPIAUCSL", label = "CPI: All Urban Consumers"),
  list(col = "CPIAPPSL", label = "CPI: Apparel"),
  list(col = "CPITRNSL", label = "CPI: Transportation"),
  list(col = "CPIMEDSL", label = "CPI: Medical Care"),
  list(col = "CUSR0000SAC", label = "CPI: Commodities ex Food/Energy"),
  list(col = "CUSR0000SAD", label = "CPI: Durables"),
  list(col = "CUSR0000SAS", label = "CPI: Services ex Energy"),
  list(col = "CPIULFSL", label = "CPI: All Items ex Food/Energy"),
  list(col = "CUSR0000SA0L2", label = "CPI: Core Goods"),
  list(col = "CUSR0000SA0L5", label = "CPI: Core Services"),
  list(col = "PCEPI", label = "PCE Price Index"),
  list(col = "DNDGRG3M086SBEA", label = "PCE: Nondurable Goods (3m)"),
  list(col = "DSERRG3M086SBEA", label = "PCE: Services (3m)"),
  list(col = "DDURRG3M086SBEA", label = "PCE: Durable Goods (3m)"),
  list(col = "OILPRICEx", label = "WTI Spot Oil"),
  list(col = "WPSID61", label = "PPI: Crude Materials"),
  list(col = "WPSID62", label = "PPI: Intermediate Manufacturing"),
  list(col = "S&P 500", label = "S&P 500"),
  list(col = "S&P div yield", label = "S&P Dividend Yield"),
  list(col = "VIXCLSx", label = "VIX")
)

fig_dir <- "./wsj_figures"
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

macro_raw <- read.csv("./2025-10-MD.csv", stringsAsFactors = FALSE, check.names = FALSE)
macro_raw$sasdate <- suppressWarnings(as.Date(macro_raw[[1]], format = "%m/%d/%Y"))
macro_raw <- macro_raw[!is.na(macro_raw$sasdate), ]

topic_dates <- df$date
topic_months <- df$month

build_monthly_target <- function(col_name) {
  if (!col_name %in% names(macro_raw)) return(NULL)
  m_dates <- macro_raw$sasdate
  m_vals <- macro_raw[[col_name]]
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
    summarise(y = mean(y, na.rm = TRUE)) %>%
    right_join(tibble(month = y_m_list), by = "month") %>%
    arrange(factor(month, levels = y_m_list))
  y_df$y
}

metrics_all <- list()

for (tgt in targets) {
  col_name <- tgt$col
  label <- tgt$label
  tag <- gsub("[^A-Za-z0-9]+", "_", tolower(col_name))


  y_raw <- build_monthly_target(col_name)
  if (is.null(y_raw)) {
    message("Skipping ", col_name, " (not found in macro file)")
    next
  }
  if (all(is.na(y_raw))) {
    message("Skipping ", col_name, " (all NA after alignment)")
    next
  }
  y_ratio <- make_ratio(y_raw)
  y <- scale(y_ratio)


  unem_df <- df %>%
    group_by(month) %>%
    summarise(y = mean(UNRATE, na.rm = TRUE))
  unem <- scale(unem_df$y)

  h <- 15
  test_idx <- (length(y) - h + 1):length(y)
  obs_idx <- setdiff(seq_along(y), test_idx)


  select_res <- select_topics(y, X_m, train_idx = obs_idx, penalty_w = 0.1)
  strong_idx <- select_res$selected
  if (length(strong_idx) == 0) {
    select_res <- select_topics(y, X_m, train_idx = obs_idx, penalty_w = 0.05)
    strong_idx <- select_res$selected
  }
  if (length(strong_idx) == 0) {
    select_res <- select_topics(y, X_m, train_idx = obs_idx, penalty_w = 0.01)
    strong_idx <- select_res$selected
  }


  ar_fit0 <- auto.arima(y[obs_idx], max.q = 0, D = 0, seasonal = FALSE,
                        allowmean = FALSE, allowdrift = FALSE)
  p_use <- max(1, ar_fit0$arma[1])
  ar_fit <- Arima(y[obs_idx], order = c(p_use, 0, 0), include.mean = FALSE)
  ar_predictions_fit <- forecast(ar_fit, h = length(test_idx))
  ar_prediction <- ar_predictions_fit$mean
  ar_interval <- cbind(ar_predictions_fit$lower[, 2], ar_predictions_fit$upper[, 2])


  ar_unem_fit <- Arima(y[obs_idx], order = c(p_use, 0, 0),
                       xreg = unem[obs_idx], include.mean = FALSE)
  ar_unem_predictions_fit <- forecast(ar_unem_fit, h = length(test_idx), xreg = unem[test_idx])
  ar_unem_predictions <- ar_unem_predictions_fit$mean
  ar_unem_interval <- cbind(ar_unem_predictions_fit$lower[, 2], ar_unem_predictions_fit$upper[, 2])


  p1 <- p_use
  p2 <- 1
  powered_lda_interval_fit <- LLM_TS.BOOT(
    X_m, y, tilde_y_m, p1, p2, obs_idx, test_idx, strong_idx, h, B = 200
  )
  LLM_prediction <- powered_lda_interval_fit[[1]]
  LLM_interval <- powered_lda_interval_fit$interval


  y_test <- y[test_idx]
  mse <- function(pred) mean((pred - y_test)^2, na.rm = TRUE)
  sign_err <- function(pred) mean(sign(pred) != sign(y_test), na.rm = TRUE)
  metrics_all[[length(metrics_all) + 1]] <- data.frame(
    target = label,
    model = c("AR", "ARX", "LLM-TS"),
    mse = c(mse(ar_prediction), mse(ar_unem_predictions), mse(LLM_prediction)),
    sign_error_rate = c(sign_err(ar_prediction), sign_err(ar_unem_predictions), sign_err(LLM_prediction))
  )


  dates <- as.Date(paste0(y_m_list, "-01"))
  forecast_start_idx <- min(test_idx)
  forecast_start_date <- dates[forecast_start_idx]

  ts_data <- data.frame(
    Date = dates,
    Actual = y
  ) %>%
    mutate(
      AR_fit = ifelse(Date >= forecast_start_date, ar_prediction, NA),
      AR_lwr = ifelse(Date >= forecast_start_date, ar_interval[, 1], NA),
      AR_upr = ifelse(Date >= forecast_start_date, ar_interval[, 2], NA),
      ARX_fit = ifelse(Date >= forecast_start_date, ar_unem_predictions, NA),
      ARX_lwr = ifelse(Date >= forecast_start_date, ar_unem_interval[, 1], NA),
      ARX_upr = ifelse(Date >= forecast_start_date, ar_unem_interval[, 2], NA),
      LLM_fit = ifelse(Date >= forecast_start_date, LLM_prediction, NA),
      LLM_lwr = ifelse(Date >= forecast_start_date, LLM_interval[, 1], NA),
      LLM_upr = ifelse(Date >= forecast_start_date, LLM_interval[, 2], NA)
    )

  his_plot_df <- ts_data[c("Actual", "AR_fit", "ARX_fit", "LLM_fit")]
  colnames(his_plot_df) <- c("Actual", "AR", "ARX", "LLM-TS")

  hist_data <- his_plot_df %>%
    pivot_longer(
      cols = everything(),
      names_to = "Type",
      values_to = "Value"
    ) %>%
    mutate(Type = factor(Type, levels = c("Actual", "AR", "ARX", "LLM-TS")))

  p_time <- ggplot(ts_data, aes(x = Date)) +
    geom_line(aes(y = Actual, color = "Actual"), linewidth = 0.8) +
    geom_line(aes(y = AR_fit, color = "AR"), linetype = "dashed", linewidth = 0.6) +
    geom_ribbon(aes(ymin = AR_lwr, ymax = AR_upr, fill = "AR"), alpha = 0.25) +
    geom_line(aes(y = ARX_fit, color = "ARX"), linewidth = 0.6) +
    geom_ribbon(aes(ymin = ARX_lwr, ymax = ARX_upr, fill = "ARX"), alpha = 0.25) +
    geom_line(aes(y = LLM_fit, color = "LLM-TS"), linewidth = 1.0) +
    geom_ribbon(aes(ymin = LLM_lwr, ymax = LLM_upr, fill = "LLM-TS"), alpha = 0.5) +
    geom_vline(xintercept = as.numeric(ts_data$Date[forecast_start_idx]),
               linetype = "dotted", color = "gray50") +
    annotate("text", x = ts_data$Date[forecast_start_idx], y = min(ts_data$Actual, na.rm = TRUE),
             label = "Forecast Start", hjust = 1.1, vjust = -0.5, size = 3, color = "gray30") +
    scale_color_manual(
      name = NULL,
      values = c(
        "Actual" = "#4D4D4D",
        "AR" = "#00BFC4",
        "ARX" = "#00A5FF",
        "LLM-TS" = "#F8766D"
      ),
      breaks = c("LLM-TS", "ARX", "AR", "Actual")
    ) +
    scale_fill_manual(
      values = c(
        "AR" = scales::alpha("#00BFC4", 0.15),
        "ARX" = scales::alpha("#00A5FF", 0.15),
        "LLM-TS" = scales::alpha("#F8766D", 0.25)
      ),
      guide = "none"
    ) +
    labs(y = label, x = NULL) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.box.margin = margin(t = -10, unit = "pt")
    )

  p_density <- ggplot(hist_data, aes(x = Value, color = Type, fill = Type, linetype = Type)) +
    geom_density(
      data = ~ filter(.x, Type == "Actual"),
      alpha = 0.3,
      linewidth = 0.8,
      bw = 0.3
    ) +
    geom_density(
      data = ~ filter(.x, Type != "Actual"),
      alpha = 0.1,
      linewidth = 0.8,
      bw = 0.3
    ) +
    geom_segment(
      aes(x = -Inf, xend = Inf, y = 0, yend = 0),
      color = "grey50",
      linewidth = 0.3
    ) +
    scale_fill_manual(
      values = c(
        "Actual" = "#4D4D4D",
        "AR" = "#00BFC4",
        "ARX" = "#00A5FF",
        "LLM-TS" = "#F8766D"
      )
    ) +
    scale_color_manual(
      values = c(
        "Actual" = "#4D4D4D",
        "AR" = "#00BFC4",
        "ARX" = "#00A5FF",
        "LLM-TS" = "#F8766D"
      )
    ) +
    scale_linetype_manual(
      values = c(
        "Actual" = "solid",
        "AR" = "dashed",
        "ARX" = "solid",
        "LLM-TS" = "solid"
      )
    ) +
    labs(x = label, y = "Density") +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 5, unit = "mm"),
      axis.line.x = element_blank()
    ) +
    guides(fill = "none", color = "none", linetype = "none")

  final_plot <- p_time / p_density +
    plot_layout(
      heights = c(1, 1),
      guides = "collect"
    ) +
    plot_annotation(theme = theme(
      legend.position = "bottom",
      legend.justification = "center",
      legend.box.margin = margin(t = -15, b = 0, unit = "pt"),
      plot.margin = margin(b = 15, unit = "pt")
    ))


  ggsave(file.path(fig_dir, paste0("wsj_plot_", tag, "_with_density.pdf")),
         final_plot, width = 10, height = 8)
  ggsave(file.path(fig_dir, paste0("wsj_plot_", tag, "_no_density.pdf")),
         p_time, width = 10, height = 5)
}

if (length(metrics_all)) {
  metrics_df <- do.call(rbind, metrics_all)
  write.csv(metrics_df, "wsj_plot_metrics.csv", row.names = FALSE)
  print(metrics_df)
}
