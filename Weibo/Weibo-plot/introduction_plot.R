setwd("/Users/sunao/Desktop/LLM-TS-Revision/code/real_data/Weibo-final/")
source('./utility.R')

library(dplyr)
library(tidyr)
library(forecast)
library(MTS)
library(lubridate)
library(ggplot2)
df <- read.csv('./lda_final.csv')
colnames(df)
df_new <- df
df_new$cpi <- df$cpi_lastmonth_100
df_new$time <- as.Date(df_new$time)
df_new$pred_scoreInflation <- df_new$aveScore
df_new <- df_new %>%
  arrange(time)

df_new <- df_new %>%
  mutate(
    month = format(time,"%Y-%m"),
    day = day(time),
    period = case_when(
      day <= 10 ~ "up",
      day <= 20 ~ "middle",
      TRUE ~ "down"
    )
  )


df_new <- df_new[,!names(df_new) %in% c("time")]

df_summary <- df_new %>%
  group_by(month, period) %>%
  summarise(across(everything(), mean, na.rm = TRUE))




y_m_list <- sort(unique(df_summary$month))
tilde_y_m <- matrix(ncol = 3, nrow = length(y_m_list))
for(i in 1:length(y_m_list)){
  res_ym <- y_m_list[i]
  tilde_y_m[i,] <- df_summary[df_summary$month == res_ym,]$pred_scoreInflation
}



topic_cols <- names(df)[grepl("^X\\d+$|^\\d+$", names(df))]


X_m <- matrix(ncol = 20, nrow = length(y_m_list))
for(i in 1:length(y_m_list)){
  res_ym <- y_m_list[i]
  X_m[i,] <- apply(df_summary[df_summary$month == res_ym,topic_cols],2,mean)
}

X_m <- scale(X_m)





y_df <- df_new %>%
  group_by(month) %>%
  summarise(y = mean(cpi, na.rm = TRUE))
y = y_df$y
y = y-100
y = scale(y)

y_df <- df_new %>%
  group_by(month) %>%
  summarise(y = mean(unem, na.rm = TRUE))
unem = y_df$y
unem = scale(unem)
h = 15
test_idx <- (length(y)-h+1):length(y)
obs_idx <- setdiff(1:length(y), test_idx)
strong_idx <- c(5,10)







ar_fit <- auto.arima(y[obs_idx],max.q = 0, D = 0,seasonal = FALSE, allowmean = FALSE, allowdrift = FALSE)
ar_fit <- Arima(y[obs_idx], order = c(ar_fit$arma[1],0,0),include.mean = FALSE)
ar_predictions_fit <- forecast(ar_fit, h = length(test_idx))
ar_prediction <- ar_predictions_fit$mean
ar_interval <- cbind(ar_predictions_fit$lower[,2], ar_predictions_fit$upper[,2])


ar_unem_fit <- Arima(y[obs_idx], order = c(ar_fit$arma[1],0,0), xreg = unem[obs_idx], include.mean = FALSE)
ar_unem_predictions_fit <- forecast(ar_unem_fit, h = length(test_idx), xreg = unem[test_idx])
ar_unem_predictions <- ar_unem_predictions_fit$mean
ar_unem_interval <- cbind(ar_unem_predictions_fit$lower[,2], ar_unem_predictions_fit$upper[,2])


p1=ar_fit$arma[1]
p2=1

powered_lda_interval_fit <- LLM_TS.BJ(X_m, y, tilde_y_m, p1, p2, obs_idx, test_idx, 
                                      strong_idx, h)

LLM_prediction <- powered_lda_interval_fit[[1]]
LLM_interval <- powered_lda_interval_fit$interval


dates <- seq.Date(from = as.Date("2019-01-01"), 
                  to = as.Date("2025-09-01"), 
                  by = "month")
n <- length(dates)
actual <- y


ts_data <- data.frame(
  Date = dates,
  Actual = actual
)

forecast_start_date <- dates[max(obs_idx)+1]



ts_data <- ts_data %>%
  mutate(

    AR_fit = ifelse(Date >= forecast_start_date, ar_prediction, NA),
    AR_lwr = ifelse(Date >= forecast_start_date, ar_interval[,1], NA),
    AR_upr = ifelse(Date >= forecast_start_date, ar_interval[,2], NA),
    

    ARX_fit = ifelse(Date >= forecast_start_date, ar_unem_predictions, NA),
    ARX_lwr = ifelse(Date >= forecast_start_date, ar_unem_interval[,1], NA),
    ARX_upr = ifelse(Date >= forecast_start_date, ar_unem_interval[,2], NA),
    

    LLM_fit = ifelse(Date >= forecast_start_date, LLM_prediction, NA),
    LLM_lwr = ifelse(Date >= forecast_start_date, LLM_interval[,1], NA),
    LLM_upr = ifelse(Date >= forecast_start_date, LLM_interval[,2], NA)
  )




his_plot_df <- ts_data[46:nrow(ts_data), c('Actual', 'AR_fit', 'ARX_fit', 'LLM_fit')]
colnames(his_plot_df) <- c("Actual", "AR", "ARX", "LLM-CPI")
long_data <- pivot_longer(his_plot_df, cols = c("Actual", "AR", "ARX", "LLM-CPI"), 
                          names_to = "Type", values_to = "Value")







library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)






hist_data <- his_plot_df %>%
  pivot_longer(
    cols = everything(),
    names_to = "Type",
    values_to = "Value"
  ) %>%
  mutate(Type = factor(Type, levels = c("Actual", "AR", "ARX", "LLM-CPI")))





p_time <- ggplot(ts_data, aes(x = Date)) +
  geom_line(aes(y = Actual, color = "Actual"), linewidth = 0.8) +
  geom_line(aes(y = AR_fit, color = "AR"), linetype = "dashed", linewidth = 0.6) +
  geom_ribbon(aes(ymin = AR_lwr, ymax = AR_upr, fill = "AR"), alpha = 0.25) +
  geom_line(aes(y = ARX_fit, color = "ARX"), linewidth = 0.6) +
  geom_ribbon(aes(ymin = ARX_lwr, ymax = AR_upr, fill = "ARX"), alpha = 0.25) +
  geom_line(aes(y = LLM_fit, color = "LLM-CPI"), linewidth = 1.0) +
  geom_ribbon(aes(ymin = LLM_lwr, ymax = LLM_upr, fill = "LLM-CPI"), alpha = 0.5) +
  geom_vline(xintercept = as.numeric(ts_data$Date[46]), 
             linetype = "dotted", color = "gray50") +
  annotate("text", x = ts_data$Date[46], y = min(ts_data$Actual, na.rm = TRUE),
           label = "Forecast Start", hjust = 1.1, vjust = -0.5, size = 3, color = "gray30") +
  scale_color_manual(
    name = NULL,
    values = c(
      "Actual" = "#4D4D4D",
      "AR" = "#00BFC4",
      "ARX" = "#00A5FF",
      "LLM-CPI" = "#F8766D"
    ),
    breaks = c("LLM-CPI", "ARX", "AR", "Actual")
  ) +
  scale_fill_manual(
    values = c(
      "AR" = scales::alpha("#00BFC4", 0.15),
      "ARX" = scales::alpha("#00A5FF", 0.15),
      "LLM-CPI" = scales::alpha("#F8766D", 0.25)
    ),
    guide = "none"
  ) +
  labs(y = "CPI", x = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.box.margin = margin(t = -10, unit = "pt")
  )





p_density <- ggplot(hist_data, aes(x = Value, color = Type, fill = Type)) +
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
      "LLM-CPI" = "#F8766D"
    )
  ) +
  scale_color_manual(
    values = c(
      "Actual" = "#4D4D4D",
      "AR" = "#00BFC4",
      "ARX" = "#00A5FF",
      "LLM-CPI" = "#F8766D"
    )
  ) +
  labs(x = "CPI", y = "Density") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 5, unit = "mm"),
    axis.line.x = element_blank()
  ) +
  guides(fill = "none", color = "none")





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


final_plot
