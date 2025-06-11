rm(list = ls())
gc()
source('./utility.R')
df <- read.csv('../meanTopicDistri_mergeCPI100_Unemployment_df19-23_20250301.csv')
colnames(df)


df_new <- df[,c(1,5:24)] 
df_new$time <- as.Date(df_new$time) # transfer variable ''time'' to Date type
df_new <- df_new %>%
  arrange(time) 

df_new <- df_new %>%
  mutate(
    month = format(time,"%Y-%m"),  # extract year and month
  )


df_new <- df_new[,!names(df_new) %in% c("time")]

df_summary <- df_new %>%
  group_by(month) %>%
  summarise(across(everything(), mean, na.rm = TRUE))


X_m <- df_summary[,2:ncol(df_summary)]
X_m <- as.matrix(X_m)
X_m <- scale(X_m)
X_select <- X_m[,c(11,12)]

######## DATA GENERATING PROCESS #############

phi1 <- c(0.5, -0.3)  # ARX(2) lag coefficients
#phi2 <- c(0.2, 0.2)
phi2 <- matrix(c(rep(0.2, 3) ,rep(-0.2,3), rep(-0.1,3)), ncol = 3, byrow = TRUE) #VARX(1) lag coefficients
beta1 <- c(0.7, -0.2)  #exogenous variable for ARX(2)

beta2 <- matrix(c(rep(0.1, 3),rep(-0.1,3)), ncol = 2, byrow = TRUE) #exogenous variable for VARX

generate_arx_data <- function(X, rho) {
  #n: length of time series
  #phi1: the AR coefficients \in \R^{q_1} for target
  #beta1: the exogenous coefficients for target
  #phi2: the VAR coefficients matrix \in \R^{q_2 \times q_2} for surrogate model
  #beta2: the exogenous coefficients matrix for surrogate model
  #rho: the correlations between random errors of target and surrogates.
  n <- nrow(X)
  p1 <- length(phi1)  
  p2 <- 1
  k <- length(beta1) # dimension of exogenous variable
  r <- ncol(phi2)  #dimension of surrogate vector
  
  sigma = rho * t(matrix(rep(1,r+1), nrow = 1)) %*% matrix(rep(1,r+1), nrow = 1)
  diag(sigma) = 1
  #sigma <- matrix(c(1, rho, rho, 1), nrow = 2) 
  error_c <- 0.1*mvrnorm(n, mu = rep(0,r+1), Sigma = sigma) 
  
  
  y <- numeric(n)
  # 生成ARX数据
  for (t in (p1 + 1):n) {
    # AR部分
    ar_term <- sum(phi1 * y[(t - 1):(t - p1)])
    
    # X部分
    x_term <- sum(beta1 * X[t, ])
    
    # 生成y_t
    y[t] <- ar_term + x_term + error_c[t,1]
  }
  
  tilde_y <- matrix(0, nrow = n, ncol = r)
  # 生成ARX数据
  for (t in 2:n) {
    # AR部分
    ar_term <- phi2 %*% tilde_y[t-1,]
    
    # X部分
    x_term <- beta2 %*% t(X[t,,drop=FALSE])
    
    # 生成y_t
    tilde_y[t,] <- ar_term + x_term + error_c[t,2:ncol(error_c)]
  }
  
  # 返回结果
  return(list(y = y, tilde_y=tilde_y))
}

########################### PREDICTION SIMULATION ######################
h_c <- seq(8,15,1)
#h_c <- seq(2,9,1)
strong_idx <- c(11,12,1) #know the true index
p1=2
p2=1
#rho_c <- seq(0,0.9, length.out = 4)
rho_c <- seq(0.1,0.4, length.out = 4)
simu_time <- 1000
MSE_result <- matrix(nrow = 0, ncol = length(h_c))
Sign_result <- matrix(nrow = 0, ncol = length(h_c))

set.seed(1234)
for(d in 1:length(rho_c)){
  
  Res_m <- matrix(0,nrow = length(h_c), ncol = 8)
  for(simu in 1:simu_time){
    ## generate data
    
    simu_data <- generate_arx_data(X_select,rho=rho_c[d])
    y <- simu_data$y
    tilde_y_m <- simu_data$tilde_y
    
    ## h step prediction
    for(k in 1:length(h_c)){
      h <- h_c[k]
      res_c <- NULL
      test_idx <- (length(y)-h+1):length(y)
      obs_idx <- setdiff(1:length(y), test_idx)
      
      
      # AR
      #ar_fit <- auto.arima(y[obs_idx],max.q = 0, D = 0,seasonal = FALSE, allowmean = FALSE, allowdrift = FALSE)
      ar_fit <- arima(y[obs_idx], order = c(p1,0,0), include.mean = FALSE)
      ar_prediction <- predict(ar_fit, n.ahead = h)$pred
      ar_mse <- mean((ar_prediction-y[test_idx])^2)
      ar_sign <- mean(sign(ar_prediction) !=sign(y[test_idx]))
      res_c <- c(res_c, mean((ar_prediction-y[test_idx])^2))
      res_c <- c(res_c, mean(sign(ar_prediction) !=sign(y[test_idx])))
      
      #random walk
      mean_prediction <- rep(y[max(obs_idx)], length(test_idx))
      res_c <- c(res_c, mean((mean_prediction-y[test_idx])^2))
      res_c <- c(res_c, mean(sign(mean_prediction) !=sign(y[test_idx])))
      
      #Average 
      ave_prediction <- Average_prediction(y[obs_idx], length(test_idx))
      res_c <- c(res_c, mean((ave_prediction-y[test_idx])^2))
      res_c <- c(res_c, mean(sign(ave_prediction) !=sign(y[test_idx])))
      
      
      # LLM powered: lag term + LDA embedding
      powered_prediction <- Prediction_powered_ts(X_m, y, tilde_y_m, p1, p2, obs_idx, test_idx, strong_idx, h)
      
      res_c <- c(res_c, mean((powered_prediction-y[test_idx])^2))
      res_c <- c(res_c, mean(sign(powered_prediction) !=sign(y[test_idx])))
      
      Res_m[k,] <- Res_m[k,] + res_c
    }
  }
  
  Res_m <- Res_m/simu_time
  MSE_result <- rbind(MSE_result, t(Res_m[,seq(1,ncol(Res_m),by=2)]))
  Sign_result <- rbind(Sign_result,t(Res_m[,seq(2,ncol(Res_m),by=2)]))
}

rownames(MSE_result) <- rep(rho_c,each = 4) # 4 is the number of methods
colnames(MSE_result) <- h_c
rownames(Sign_result) <- rep(rho_c,each = 4)
colnames(Sign_result) <- h_c

write.csv(MSE_result, './result/over_MSE_result.csv', row.names = FALSE)
write.csv(Sign_result, './result/over_Sign_result.csv', row.names = FALSE)
MSE_result
Sign_result

