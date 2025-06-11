##new dataset
library(dplyr)
library(forecast)
library(MTS)
library(lubridate)
library(ggplot2)
library(MASS)
setwd("/Users/sunao/Desktop/prediction powered time series/code")
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



X_m <- matrix(ncol = 20, nrow = length(y_m_list))
for(i in 1:length(y_m_list)){
  res_ym <- y_m_list[i]
  X_m[i,] <- apply(df_summary[df_summary$month == res_ym,4:23],2,mean)
}


y_df <- df_new %>%
  group_by(month) %>%
  summarise(y = mean(cpi))
y = y_df$y
y = y-100
y = scale(y)


h <- 5
test_idx <- (length(y)-h+1):length(y)
obs_idx <- setdiff(1:length(y), test_idx)
ar_fit <- auto.arima(y[obs_idx],max.q = 0, D = 0,seasonal = FALSE, allowmean = FALSE, allowdrift = FALSE)
ar_fit <- arima(y[obs_idx], order = c(ar_fit$arma[1],0,0), include.mean = FALSE)
strong_idx <- order(abs(cov(ar_fit$residuals, X_m[obs_idx,])), decreasing = TRUE)[1:2]
arx_fit <- arima(y[obs_idx], order = c(2,0,0), xreg = X_m[obs_idx,strong_idx], include.mean = FALSE)
arx_predictions <- predict(arx_fit, n.ahead = h, newxreg = X_m[test_idx,strong_idx])$pred
tilde_y_fit <- VARX(zt = tilde_y_m, p=1, xt = X_m[,strong_idx], m=0, include.mean = FALSE)

tilde_y_residual_m <- tilde_y_fit$residuals[obs_idx,]

plot(density(arx_fit$residuals))
plot(density(tilde_y_residual_m[,1]))
plot(density(tilde_y_residual_m[,2]))
plot(density(tilde_y_residual_m[,3]))


library(ggplot2)
library(dplyr)

# Generate example data
set.seed(123)
data1 <- data.frame(x = as.numeric(arx_fit$residuals), y = as.numeric(tilde_y_residual_m[,1]), group = "First")
data2 <- data.frame(x = as.numeric(arx_fit$residuals), y = as.numeric(tilde_y_residual_m[,2]), group = "Middle")
data3 <- data.frame(x = as.numeric(arx_fit$residuals), y = as.numeric(tilde_y_residual_m[,3]), group = "Last")

# Combine data
combined_data <- bind_rows(data1, data2, data3)

# Convert 'group' to a factor with specified order
combined_data$group <- factor(combined_data$group, levels = c("First", "Middle", "Last"))

# Remove outliers (your existing function)
remove_outliers <- function(data, x_col, y_col) {
  x_Q1 <- quantile(data[[x_col]], 0.25)
  x_Q3 <- quantile(data[[x_col]], 0.75)
  x_IQR <- x_Q3 - x_Q1
  x_lower_bound <- x_Q1 - 1.5 * x_IQR
  x_upper_bound <- x_Q3 + 1.5 * x_IQR
  
  y_Q1 <- quantile(data[[y_col]], 0.25)
  y_Q3 <- quantile(data[[y_col]], 0.75)
  y_IQR <- y_Q3 - y_Q1
  y_lower_bound <- y_Q1 - 1.5 * y_IQR
  y_upper_bound <- y_Q3 + 1.5 * y_IQR
  
  data_cleaned <- data %>%
    filter(data[[x_col]] >= x_lower_bound & data[[x_col]] <= x_upper_bound &
             data[[y_col]] >= y_lower_bound & data[[y_col]] <= y_upper_bound)
  
  return(data_cleaned)
}

# Apply outlier removal
combined_data_cleaned <- combined_data %>%
  group_by(group) %>%
  group_modify(~ remove_outliers(.x, "x", "y"))

# Plot with facets in the desired order
ggplot(combined_data_cleaned, aes(x = x, y = y)) +
  geom_density_2d(aes(color = ..level..), size = 1) +
  geom_point(alpha = 0.3, color = "blue") +
  scale_color_viridis_c() +
  facet_wrap(~ group, ncol = 3) +  # Will now follow factor levels
  labs(x = "LLM inflation index", y = "CPI", color = "Density level") +
  xlim(-3, 3) + ylim(-0.075, 0.075) +
  theme_minimal()

