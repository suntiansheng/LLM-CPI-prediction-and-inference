setwd("/Users/sunao/Desktop/prediction powered time series/code")
source('./final code/utility.R')
##new dataset
library(dplyr)
library(forecast)
library(MTS)
library(lubridate)
library(ggplot2)
df <- read.csv('./meanTopicDistri_mergeCPI100_Unemployment_df19-23_20250301.csv')
colnames(df)
df_new <- df[,c(1,3,5:24)]
df_new$unem <-df$全国城镇调查失业率...
df_new$cpi <- df$cpi_lastmonth_100
df_new$time <- as.Date(df_new$time)
df_new <- df_new %>%
  arrange(time)

df_new <- df_new %>%
  mutate(
    month = format(time,"%Y-%m"),  # 提取年月
    day = day(time),                # 提取日
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


#max_num_day <- min(table(df_new$year_month))

y_m_list <- sort(unique(df_summary$month))
tilde_y_m <- matrix(ncol = 3, nrow = length(y_m_list))
for(i in 1:length(y_m_list)){
  res_ym <- y_m_list[i]
  tilde_y_m[i,] <- df_summary[df_summary$month == res_ym,]$pred_scoreInflation
}

#tilde_y_m <- scale(tilde_y_m-0.5, center = FALSE)



X_m <- matrix(ncol = 20, nrow = length(y_m_list))
for(i in 1:length(y_m_list)){
  res_ym <- y_m_list[i]
  X_m[i,] <- apply(df_summary[df_summary$month == res_ym,4:23],2,mean)
}

#X_m <- scale(X_m, scale = FALSE)


#X_lda <- X_m[,c(11,12)]
#write.csv(X_lda,file = 'lda_embedding.csv', row.names = FALSE)

y_df <- df_new %>%
  group_by(month) %>%
  summarise(y = mean(cpi))
y = y_df$y
y = y-100
y = scale(y)

y_df <- df_new %>%
  group_by(month) %>%
  summarise(y = mean(unem))
unem = y_df$y
unem = scale(unem)
h = 15
test_idx <- (length(y)-h+1):length(y)
obs_idx <- setdiff(1:length(y), test_idx)
strong_idx <- c(11,12)


#### Inference analysis



#AR
ar_fit <- auto.arima(y[obs_idx],max.q = 0, D = 0,seasonal = FALSE, allowmean = FALSE, allowdrift = FALSE)
ar_fit <- Arima(y[obs_idx], order = c(ar_fit$arma[1],0,0),include.mean = FALSE)
ar_predictions_fit <- forecast(ar_fit, h = length(test_idx))
ar_prediction <- ar_predictions_fit$mean
ar_interval <- cbind(ar_predictions_fit$lower[,2], ar_predictions_fit$upper[,2])

# ARX
ar_unem_fit <- Arima(y[obs_idx], order = c(ar_fit$arma[1],0,0), xreg = unem[obs_idx], include.mean = FALSE)
ar_unem_predictions_fit <- forecast(ar_unem_fit, h = length(test_idx), xreg = unem[test_idx])
ar_unem_predictions <- ar_unem_predictions_fit$mean
ar_unem_interval <- cbind(ar_unem_predictions_fit$lower[,2], ar_unem_predictions_fit$upper[,2])

# LLM-CPI
p1=ar_fit$arma[1]
p2=1

powered_lda_interval_fit <- Prediction_BJ(X_m, y, tilde_y_m, p1, p2, obs_idx, test_idx, 
                                      strong_idx, h)

LLM_prediction <- powered_lda_interval_fit[[1]]
LLM_interval <- powered_lda_interval_fit$interval


############### Density plot #####################
# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)  # For combining plots

# Time Series Plot
p1 <- ggplot(ts_data, aes(x = Date)) +
  # Actual data line
  geom_line(aes(y = Actual, color = "Actual"), linewidth = 1) +
  
  # AR fit and ribbon
  geom_line(aes(y = AR_fit, color = "AR"), linetype = "dashed", linewidth = 0.7, na.rm = TRUE) +
  geom_ribbon(aes(ymin = AR_lwr, ymax = AR_upr, fill = "AR"), alpha = 0.2, na.rm = TRUE) +
  
  # ARX fit and ribbon
  geom_line(aes(y = ARX_fit, color = "ARX"), linewidth = 0.7, na.rm = TRUE) +
  geom_ribbon(aes(ymin = ARX_lwr, ymax = ARX_upr, fill = "ARX"), alpha = 0.2, na.rm = TRUE) +
  
  # LLM fit and ribbon
  geom_line(aes(y = LLM_fit, color = "LLM-CPI"), linewidth = 1.5, na.rm = TRUE) +
  geom_ribbon(aes(ymin = LLM_lwr, ymax = LLM_upr, fill = "LLM-CPI"), alpha = 0.4, color = "#e31a1c", linewidth = 0.3, na.rm = TRUE) +
  
  # Forecast start marker
  geom_vline(xintercept = forecast_start_date, linetype = "dashed", color = "gray50") +
  annotate("text", x = forecast_start_date, y = min(ts_data$Actual, na.rm = TRUE), 
           label = "Forecast start", hjust = 1.1, vjust = -1, color = "gray50") +
  
  # Custom color scheme
  scale_color_manual(
    name = "Forecasts",
    values = c("Actual" = "darkgrey", "AR" = "cyan", "ARX" = "deepskyblue", "LLM-CPI" = "red"),
    breaks = c("LLM-CPI", "ARX", "AR", "Actual")
  ) +
  scale_fill_manual(
    name = "Prediction Intervals",
    values = c("AR" = "cyan", "ARX" = "deepskyblue", "LLM-CPI" = "#e31a1c"),
    breaks = c("LLM-CPI", "ARX", "AR")
  ) +
  labs(
    y = "CPI",
    x = "Date"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    legend.title = element_text(face = "bold")
  )

# Histogram Plot
his_plot_df <- ts_data[46:nrow(ts_data), c('Actual', 'AR_fit', 'ARX_fit', 'LLM_fit')]
colnames(his_plot_df) <- c("Actual", "AR", "ARX", "LLM-CPI")
long_data <- pivot_longer(his_plot_df, cols = c("Actual", "AR", "ARX", "LLM-CPI"), 
                          names_to = "Type", values_to = "Value")

bin_num <- 7

p2 <- ggplot(long_data, aes(x = Value, fill = Type)) +
  # Add histograms for each type with transparency
  geom_histogram(
    data = filter(long_data, Type == "Actual"),
    bins = bin_num,
    alpha = 0.6,
    position = "identity",
    color = "white"
  ) +
  geom_histogram(
    data = filter(long_data, Type == "AR"),
    bins = bin_num,
    alpha = 0.2,
    position = "identity",
    color = "white"
  ) +
  geom_histogram(
    data = filter(long_data, Type == "ARX"),
    bins = bin_num,
    alpha = 0.2,
    position = "identity",
    color = "white"
  ) +
  geom_histogram(
    data = filter(long_data, Type == "LLM-CPI"),
    bins = bin_num,
    alpha = 0.6,
    position = "identity",
    color = "white"
  ) +
  # Define custom colors for each type and reorder the legend
  scale_fill_manual(
    values = c(
      "LLM-CPI" = "red",
      "ARX" = "deepskyblue",
      "AR" = "cyan",
      "Actual" = "darkgrey"
    ),
    name = "Type",
    breaks = c("LLM-CPI", "ARX", "AR", "Actual")  # Reorder the legend
  ) +
  labs(
    x = "CPI",
    y = NULL
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    legend.title = element_text(face = "bold")
  )

# Combine the two plots in a single column
combined_plot <- p1 / p2

# Show the combined plot
combined_plot


##############
# Load Required Packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# ------------------------------------------
# Step 1: Prepare Data
# ------------------------------------------

# Convert to long format for histogram
hist_data <- his_plot_df %>%
  pivot_longer(
    cols = everything(),
    names_to = "Type",
    values_to = "Value"
  ) %>%
  mutate(Type = factor(Type, levels = c("Actual", "AR", "ARX", "LLM-CPI")))

# ------------------------------------------
# Step 2: Create Time Series Plot without Legend Title
# ------------------------------------------

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
    name = NULL,  # Removed legend title
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

# ------------------------------------------
# Step 3: Create Histogram Plot with X-axis Line
# ------------------------------------------

p_hist <- ggplot(hist_data, aes(x = Value)) +
  geom_histogram(
    data = ~ filter(.x, Type == "Actual"),
    aes(fill = Type),
    bins = 7,
    alpha = 0.6,
    color = "white"
  ) +
  geom_histogram(
    data = ~ filter(.x, Type != "Actual"),
    aes(color = Type),
    bins = 7,
    fill = NA,
    linewidth = 0.6
  ) +
  # Add grey x-axis line
  geom_segment(
    aes(x = -Inf, xend = Inf, y = 0, yend = 0),
    color = "grey50",
    linewidth = 0.3
  ) +
  scale_fill_manual(values = c("Actual" = "#4D4D4D")) +
  scale_color_manual(
    values = c(
      "AR" = "#00BFC4",
      "ARX" = "#00A5FF",
      "LLM-CPI" = "#F8766D"
    )
  ) +
  labs(x = "CPI Values", y = "Count") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 5, unit = "mm"),
    axis.line.x = element_blank()  # Remove default axis line
  ) +
  guides(fill = "none", color = "none")

# ------------------------------------------
# Step 4: Combine Plots
# ------------------------------------------

final_plot <- p_time / p_hist +
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

# Display final plot
final_plot


############# KDE
# Load Required Packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# ------------------------------------------
# Step 1: Prepare Data
# ------------------------------------------

# Convert to long format for density plot
hist_data <- his_plot_df %>%
  pivot_longer(
    cols = everything(),
    names_to = "Type",
    values_to = "Value"
  ) %>%
  mutate(Type = factor(Type, levels = c("Actual", "AR", "ARX", "LLM-CPI")))

# ------------------------------------------
# Step 2: Create Time Series Plot without Legend Title
# ------------------------------------------

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
    name = NULL,  # Removed legend title
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

# ------------------------------------------
# Step 3: Create Smoother Kernel Density Plot with X-axis Line
# ------------------------------------------

p_density <- ggplot(hist_data, aes(x = Value, color = Type, fill = Type)) +
  geom_density(
    data = ~ filter(.x, Type == "Actual"),
    alpha = 0.3,
    linewidth = 0.8,
    bw = 0.3  # Larger bandwidth for smoother lines
  ) +
  geom_density(
    data = ~ filter(.x, Type != "Actual"),
    alpha = 0.1,
    linewidth = 0.8,
    bw = 0.3  # Larger bandwidth for smoother lines
  ) +
  # Add grey x-axis line
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
    axis.line.x = element_blank()  # Remove default axis line
  ) +
  guides(fill = "none", color = "none")

# ------------------------------------------
# Step 4: Combine Plots
# ------------------------------------------

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

# Display final plot
final_plot
