setwd("/Users/sunao/Desktop/prediction powered time series/code/final code")
source('./utility.R')
library(dplyr)
library(forecast)
library(MTS)
library(lubridate)


############### 19-21 ##################
##new dataset
df <- read.csv('../intermediate_code/3year_meanTopic20Distri_mergeCPI_df19-21_20250314.csv')
colnames(df)
#3year_meanTopic20Distri_mergeCPI_df19-21_20250314.csv
#3year_meanTopic20Distri_mergeCPI_df21-23_20250314.csv
Top_num <- 20
colnames(df)
df_new <- df[,c(1,3,7:(Top_num+6))]
colnames(df_new)
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

tilde_y_m <- scale(tilde_y_m-0.5, center = FALSE)
#ilde_y_m <- scale(tilde_y_m)

X_m <- matrix(ncol = Top_num, nrow = length(y_m_list))
for(i in 1:length(y_m_list)){
  res_ym <- y_m_list[i]
  X_m[i,] <- apply(df_summary[df_summary$month == res_ym,4:(Top_num+3)],2,mean)
}
X_m <- scale(X_m, scale = FALSE)

#X_lda <- X_m[,c(11,12)]
#write.csv(X_lda,file = 'lda_embedding.csv', row.names = FALSE)

y_df <- df_new %>%
  group_by(month) %>%
  summarise(y = mean(cpi))
y = y_df$y
y = y-100
y = scale(y,center = FALSE)

#tilde_y_fit <- VARX(zt = tilde_y_m, p=1, xt = X_m[,strong_idx], m=0, include.mean = FALSE)
#arx_fit <- arima(y, order = c(ar_fit$arma[1],0,0), xreg = X_m[,strong_idx], include.mean = FALSE)
#cor(cbind(arx_fit$residuals[2:length(arx_fit$residuals)],tilde_y_fit$residuals))

#############################################################################
#############################################################################
set.seed(1234)
h_c <- seq(2,6,1)
res_c <- NULL
rank_c <- NULL
for(k in 1:length(h_c)){
  h <- h_c[k]
  test_idx <- (length(y)-h+1):length(y)
  obs_idx <- setdiff(1:length(y), test_idx)
  ar_fit <- auto.arima(y[obs_idx], max.q = 0, D = 0,seasonal = FALSE, allowmean = FALSE, allowdrift = FALSE)
  ar_fit <- arima(y[obs_idx], order = c(ar_fit$arma[1],0,0), include.mean = FALSE)
  ar_prediction <- predict(ar_fit, n.ahead = h)$pred
  
  mse_c <- NULL
  for(i in 1:Top_num){
    strong_idx <- order(abs(cov(ar_fit$residuals, X_m[obs_idx,])), decreasing = TRUE)[1:i]
    arx_fit <- arima(y[obs_idx], order = c(ar_fit$arma[1],0,0), xreg = X_m[obs_idx,strong_idx], include.mean = FALSE)
    arx_predictions <- predict(arx_fit, n.ahead = h, newxreg = X_m[test_idx,strong_idx])$pred
    
    mse_c[i] <- mean((arx_predictions-y[test_idx])^2)
    
  }
  rank_c[k] <- which.min(mse_c)
}
names(which.max(table(rank_c)))

num_topic <- as.numeric(names(which.max(table(rank_c))))
strong_idx <- order(abs(cov(ar_fit$residuals, X_m[obs_idx,])), decreasing = TRUE)[1:num_topic]
# The final analysis

#h_c <- seq(8,15,1)
set.seed(1234)
Res_m <- matrix(nrow = length(h_c), ncol = 11)
for(k in 1:length(h_c)){
  h <- h_c[k]
  res_c <- NULL
  test_idx <- (length(y)-h+1):length(y)
  obs_idx <- setdiff(1:length(y), test_idx)
  
  #AR
  ar_fit <- auto.arima(y[obs_idx],max.q = 0, D = 0,seasonal = FALSE, allowmean = FALSE, allowdrift = FALSE)
  ar_fit <- arima(y[obs_idx], order = c(ar_fit$arma[1],0,0), include.mean = FALSE)
  ar_prediction <- predict(ar_fit, n.ahead = h)$pred
  
  res_c <- c(res_c, mean((ar_prediction-y[test_idx])^2))
  res_c <- c(res_c, mean(sign(ar_prediction) !=sign(y[test_idx])))
  
  #random walk
  mean_prediction <- rep(y[max(obs_idx)], length(test_idx))
  res_c <- c(res_c, sqrt(mean((mean_prediction-y[test_idx])^2)))
  res_c <- c(res_c, mean(sign(mean_prediction) !=sign(y[test_idx])))
  
  #Average 
  ave_prediction <- Average_prediction(y[obs_idx], length(test_idx))
  res_c <- c(res_c, sqrt(mean((ave_prediction-y[test_idx])^2)))
  res_c <- c(res_c, mean(sign(ave_prediction) !=sign(y[test_idx])))
  
  #LDA 
  arx_fit <- arima(y[obs_idx], order = c(ar_fit$arma[1],0,0), xreg = X_m[obs_idx,strong_idx], include.mean = FALSE)
  arx_predictions <- predict(arx_fit, n.ahead = h, newxreg = X_m[test_idx,strong_idx])$pred
  
  res_c <- c(res_c, mean((arx_predictions-y[test_idx])^2))
  res_c <- c(res_c, mean(sign(arx_predictions) !=sign(y[test_idx])))
  
  # Power
  p1=ar_fit$arma[1]
  p2=1
  powered_prediction <- Prediction_powered_ts(X_m, y, tilde_y_m, p1, p2, obs_idx, test_idx, 
                                              strong_idx, h)
  
  res_c <- c(res_c, mean((powered_prediction-y[test_idx])^2))
  res_c <- c(res_c, mean(sign(powered_prediction) !=sign(y[test_idx])))
  Res_m[k,] <- c(h,res_c)
}

colnames(Res_m) <- c('h','ar_mse', 'ar_sign', 'RW_mse', 'RW_sign', 'AVE_mse','AVE_sign', 'LDA_mse', 'LDA_sign', 'LLM_LDA_mse','LLM_LDA_sign')
Res_m
#colMeans(Res_m)
#rbind(Res_m, colMeans(Res_m))
row.names(Res_m) <- Res_m[,1]
Res_m <- Res_m[,2:ncol(Res_m)]

Res_m[,seq(1,ncol(Res_m),2)]
Res_m[,seq(2,ncol(Res_m),2)]
  
#sweep(Res_m[,seq(1,ncol(Res_m),2)], 1, Res_m[,1,drop=FALSE], "/")
#sweep(Res_m[,seq(2,ncol(Res_m),2)], 1, Res_m[,2,drop=FALSE], "/")
#Res_m <- Res_m[,2:ncol(Res_m)]
#sweep(Res_m, 1, Res_m[,1,drop=FALSE], "/")
#Res_m[,seq(1,ncol(Res_m),2)]
#Res_m[,c(1,3,5,7,9,11)]
#colnames(df_summary)[strong_idx+3]

library(corrplot)
library(RColorBrewer)

col <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100)

# Create the plot with beautiful styling
corrplot(cor(cbind(tilde_y_m,y)),
         method = "circle",        # Visualization method
         type = "upper",           # Show only upper triangle
         order = "hclust",         # Cluster variables
         col = col,                # Custom color palette
         bg = "white",             # Background color
         tl.col = "darkblue",      # Text label color
         tl.cex = 0.9,             # Text label size
         tl.srt = 45,              # Text label rotation
         addCoef.col = "black",    # Correlation coefficient color
         number.cex = 0.7,         # Correlation coefficient size
         diag = FALSE,            # Hide diagonal
         outline = "white",        # Circle outline color
         mar = c(0, 0, 1, 0),      # Margins
         cl.ratio = 0.2,           # Color legend width ratio
         cl.align.text = "l",      # Color legend text alignment
         cl.offset = 0.5,          # Color legend offset
         addgrid.col = "gray90")   # Grid color

# Add title
title("Correlation Matrix of mtcars Variables", 
      line = -1, 
      cex.main = 1.2, 
      col.main = "darkblue")
#########################################################################################################
#########################################################################################################
#########################################################################################################


df <- read.csv('../intermediate_code/3year_meanTopic20Distri_mergeCPI_df21-23_20250314.csv')
colnames(df)
#3year_meanTopic20Distri_mergeCPI_df19-21_20250314.csv
#3year_meanTopic20Distri_mergeCPI_df21-23_20250314.csv
Top_num <- 20
colnames(df)
df_new <- df[,c(1,3,7:(Top_num+6))]
colnames(df_new)
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

tilde_y_m <- scale(tilde_y_m-0.5, center = FALSE)
#ilde_y_m <- scale(tilde_y_m)

X_m <- matrix(ncol = Top_num, nrow = length(y_m_list))
for(i in 1:length(y_m_list)){
  res_ym <- y_m_list[i]
  X_m[i,] <- apply(df_summary[df_summary$month == res_ym,4:(Top_num+3)],2,mean)
}
X_m <- scale(X_m, scale = FALSE)

#X_lda <- X_m[,c(11,12)]
#write.csv(X_lda,file = 'lda_embedding.csv', row.names = FALSE)

y_df <- df_new %>%
  group_by(month) %>%
  summarise(y = mean(cpi))
y = y_df$y
y = y-100
y = scale(y,center = FALSE)

#tilde_y_fit <- VARX(zt = tilde_y_m, p=1, xt = X_m[,strong_idx], m=0, include.mean = FALSE)
#arx_fit <- arima(y, order = c(ar_fit$arma[1],0,0), xreg = X_m[,strong_idx], include.mean = FALSE)
#cor(cbind(arx_fit$residuals[2:length(arx_fit$residuals)],tilde_y_fit$residuals))

#############################################################################
#############################################################################
set.seed(1234)
h_c <- seq(2,6,1)
res_c <- NULL
rank_c <- NULL
for(k in 1:length(h_c)){
  h <- h_c[k]
  test_idx <- (length(y)-h+1):length(y)
  obs_idx <- setdiff(1:length(y), test_idx)
  ar_fit <- auto.arima(y[obs_idx], max.q = 0, D = 0,seasonal = FALSE, allowmean = FALSE, allowdrift = FALSE)
  ar_fit <- arima(y[obs_idx], order = c(ar_fit$arma[1],0,0), include.mean = FALSE)
  ar_prediction <- predict(ar_fit, n.ahead = h)$pred
  
  mse_c <- NULL
  for(i in 1:Top_num){
    strong_idx <- order(abs(cor(ar_fit$residuals, X_m[obs_idx,])), decreasing = TRUE)[1:i]
    arx_fit <- arima(y[obs_idx], order = c(ar_fit$arma[1],0,0), xreg = X_m[obs_idx,strong_idx], include.mean = FALSE)
    arx_predictions <- predict(arx_fit, n.ahead = h, newxreg = X_m[test_idx,strong_idx])$pred
    
    mse_c[i] <- mean((arx_predictions-y[test_idx])^2)
    
  }
  rank_c[k] <- which.min(mse_c)
}
names(which.max(table(rank_c)))

num_topic <- as.numeric(names(which.max(table(rank_c))))
strong_idx <- order(abs(cor(ar_fit$residuals, X_m[obs_idx,])), decreasing = TRUE)[1:num_topic]

# The final analysis

#h_c <- seq(8,15,1)
set.seed(1234)
Res_m <- matrix(nrow = length(h_c), ncol = 11)
for(k in 1:length(h_c)){
  h <- h_c[k]
  res_c <- NULL
  test_idx <- (length(y)-h+1):length(y)
  obs_idx <- setdiff(1:length(y), test_idx)
  
  ar_fit <- auto.arima(y[obs_idx],max.q = 0, D = 0,seasonal = FALSE, allowmean = FALSE, allowdrift = FALSE)
  ar_fit <- arima(y[obs_idx], order = c(ar_fit$arma[1],0,0), include.mean = FALSE)
  ar_prediction <- predict(ar_fit, n.ahead = h)$pred
  
  res_c <- c(res_c, mean((ar_prediction-y[test_idx])^2))
  res_c <- c(res_c, mean(sign(ar_prediction) !=sign(y[test_idx])))
  
  #random walk
  mean_prediction <- rep(y[max(obs_idx)], length(test_idx))
  res_c <- c(res_c, sqrt(mean((mean_prediction-y[test_idx])^2)))
  res_c <- c(res_c, mean(sign(mean_prediction) !=sign(y[test_idx])))
  
  #Average 
  ave_prediction <- Average_prediction(y[obs_idx], length(test_idx))
  res_c <- c(res_c, sqrt(mean((ave_prediction-y[test_idx])^2)))
  res_c <- c(res_c, mean(sign(ave_prediction) !=sign(y[test_idx])))
  
  
  arx_fit <- arima(y[obs_idx], order = c(ar_fit$arma[1],0,0), xreg = X_m[obs_idx,strong_idx], include.mean = FALSE)
  arx_predictions <- predict(arx_fit, n.ahead = h, newxreg = X_m[test_idx,strong_idx])$pred
  
  res_c <- c(res_c, mean((arx_predictions-y[test_idx])^2))
  res_c <- c(res_c, mean(sign(arx_predictions) !=sign(y[test_idx])))
  
  
  
  
  p1=ar_fit$arma[1]
  p2=1
  powered_prediction <- Prediction_powered_ts(X_m, y, tilde_y_m, p1, p2, obs_idx, test_idx, 
                                              strong_idx, h)
  
  res_c <- c(res_c, mean((powered_prediction-y[test_idx])^2))
  res_c <- c(res_c, mean(sign(powered_prediction) !=sign(y[test_idx])))
  Res_m[k,] <- c(h,res_c)
}

colnames(Res_m) <- c('h','ar_mse', 'ar_sign', 'RW_mse', 'RW_sign', 'AVE_mse','AVE_sign', 'LDA_mse', 'LDA_sign', 'LLM_LDA_mse','LLM_LDA_sign')
Res_m
#colMeans(Res_m)
#rbind(Res_m, colMeans(Res_m))
row.names(Res_m) <- Res_m[,1]
Res_m <- Res_m[,2:ncol(Res_m)]

Res_m[,seq(1,ncol(Res_m),2)]
Res_m[,seq(2,ncol(Res_m),2)]

