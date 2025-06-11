Prediction_powered_ts <- function(X_m, y, tilde_y_m, p1, p2, obs_idx, test_idx, 
                                  strong_idx, h){
  Ti <- dim(tilde_y_m)[1]
  #tilde_y_fit <- arima(tilde_y, order = c(p2, 0,0), xreg = X)
  #tilde_y_fit <- auto.arima(tilde_y, max.q = 0, xreg = X[,strong_idx],  D = 0, seasonal = FALSE, allowmean = FALSE, allowdrift = FALSE)
  #tilde_y_fit <- VARX(zt = tilde_y_m, p=p2, xt = X_m[,strong_idx], m=0, include.mean = FALSE)
  invisible(capture.output(tilde_y_fit <- VARX(zt = tilde_y_m, p=p2, xt = X_m[,strong_idx], m=0, include.mean = FALSE)))
  #p2 <- tilde_y_fit$arma[1]
  
  #tilde_y_fit <- arima(tilde_y, order = c(p2,0,0), xreg = X[,strong_idx], include.mean = FALSE)
  
  ar_est <- tilde_y_fit$Phi
  hat_y <- matrix(0,nrow = Ti, dim(tilde_y_m)[2])
  for(t in (p2+1):Ti){
    hat_y[t,] <- (diag(ncol(tilde_y_m)) - ar_est)%*%t(tilde_y_m[t,,drop = FALSE])
  }
  
  y_obs <- y[obs_idx]
  X_aug <- cbind(X_m[,strong_idx], hat_y) # only use highly correlated index
  X_aug_obs <- X_aug[obs_idx,]
  
  #y_fit <- auto.arima(y_obs, max.q = 0, xreg = X_aug_obs, D = 0,seasonal = FALSE, allowmean = FALSE, allowdrift = FALSE)
  #p1 <- y_fit$arma[1]
  y_fit <- arima(y_obs, order = c(p1,0,0), xreg = X_aug_obs, include.mean = FALSE)
  
  #y_fit <- arima(y_obs, order = c(p1, 0,0), xreg = X_aug_obs)
  
  #future_idx <- setdiff(1:n, obs_idx)
  #X_future <- X_aug[future_idx,]
  X_future <- X_aug[test_idx,]
  #h <- length(future_idx)
  
  predictions <- predict(y_fit, n.ahead = h, newxreg = X_future)$pred
  return(predictions)
}

Prediction_BJ <- function(X_m, y, tilde_y_m, p1, p2, obs_idx, test_idx, 
                          strong_idx, h, alpha=0.05){
  #X_m \in R^{T \times p}: input exogenous predictors
  #y \in R^T: the target time series
  #tilde_y_m \in R^{T \times q}: the synthetic surrogate 
  #p1: the AR order for target ARX model
  #p2: the VAR order for synthetic surrogate 
  #obs_idx: the observed index for target data
  #test_idx: the unobserved index for target data, the union of obs_idx and test_idx is equal to the whole index
  #strong_idx: the important index for predictors
  
  
  Ti <- dim(tilde_y_m)[1]
  #tilde_y_fit <- arima(tilde_y, order = c(p2, 0,0), xreg = X)
  #tilde_y_fit <- auto.arima(tilde_y, max.q = 0, xreg = X[,strong_idx],  D = 0, seasonal = FALSE, allowmean = FALSE, allowdrift = FALSE)
  invisible(capture.output(tilde_y_fit <- VARX(zt = tilde_y_m, p=p2, xt = X_m[,strong_idx], m=0, include.mean = FALSE)))
  #p2 <- tilde_y_fit$arma[1]
  
  #tilde_y_fit <- arima(tilde_y, order = c(p2,0,0), xreg = X[,strong_idx], include.mean = FALSE)
  
  ar_est <- tilde_y_fit$Phi
  hat_y <- matrix(0,nrow = Ti, dim(tilde_y_m)[2])
  for(t in (p2+1):Ti){
    hat_y[t,] <- (diag(ncol(tilde_y_m)) - ar_est)%*%t(tilde_y_m[t,,drop = FALSE])
  }
  
  y_obs <- y[obs_idx]
  X_aug <- cbind(X_m[,strong_idx], hat_y) # only use highly correlated index
  X_aug_obs <- X_aug[obs_idx,]
  
  #y_fit <- auto.arima(y_obs, max.q = 0, xreg = X_aug_obs, D = 0,seasonal = FALSE, allowmean = FALSE, allowdrift = FALSE)
  #p1 <- y_fit$arma[1]
  y_fit <- arima(y_obs, order = c(p1,0,0), xreg = X_aug_obs, include.mean = FALSE)
  
  var_e <- mean(y_fit$residuals^2)
  alpha_hat <- y_fit$coef[1:p1]
  I_q1 <- diag(p1-1)
  I_q1 <- cbind(I_q1, 0)
  A_hat <- rbind(alpha_hat, I_q1)
  X_future <- X_aug[test_idx,]
  #h <- length(future_idx)
  predictions <- predict(y_fit, n.ahead = h, newxreg = X_future)$pred
  
  interval_m <- matrix(nrow = h, ncol = 2)
  var_h <- 0
  for(i in 1:h){
    res_var <- var_e*(A_hat%^%(i-1))[1,1]
    var_h <- res_var + var_h
    interval_m[i,1] <- predictions[i]-sqrt(var_h)*qnorm(alpha/2, lower.tail = FALSE)
    interval_m[i,2] <- predictions[i]+sqrt(var_h)*qnorm(alpha/2, lower.tail = FALSE)
  }
  return(list(predictions, interval = interval_m))
}

#library(doParallel)
#library(foreach)
#library(MTS)  # 提供VARX函数

# Bootstrap_powered_ts <- function(X_m, y, tilde_y_m, p1, p2, obs_idx, test_idx,
#                                  strong_idx, h, alpha = 0.05, B = 500) {
#   
#   # 初始化并行环境
#   cl <- makeCluster(detectCores() - 1)
#   registerDoParallel(cl)
#   on.exit(stopCluster(cl))
#   
#   # 核心计算
#   Ti <- nrow(tilde_y_m)
#   q <- ncol(tilde_y_m)
#   
#   # 拟合VARX模型
#   invisible(capture.output(
#     tilde_y_fit <- VARX(zt = tilde_y_m, p = p2, 
#                         xt = X_m[, strong_idx, drop = FALSE], 
#                         m = 0, include.mean = FALSE)
#   ))
#   
#   # 计算hat_y矩阵
#   ar_est <- tilde_y_fit$Phi
#   hat_y <- (diag(q) - ar_est) %*% t(tilde_y_m[(p2+1):Ti, ])
#   hat_y_full <- matrix(0, Ti, q)
#   hat_y_full[(p2+1):Ti, ] <- t(hat_y)
#   
#   # 准备ARIMA数据
#   X_aug <- cbind(X_m[, strong_idx, drop = FALSE], hat_y_full)
#   X_aug_obs <- X_aug[obs_idx, , drop = FALSE]
#   y_obs <- y[obs_idx]
#   
#   # 拟合ARIMA模型
#   y_fit <- arima(y_obs, order = c(p1,0,0), xreg = X_aug_obs, include.mean = FALSE)
#   
#   # 预计算预测结果（关键修正点）
#   X_future <- X_aug[test_idx, , drop = FALSE]
#   predictions <- predict(y_fit, n.ahead = h, newxreg = X_future)$pred  # 确保这行存在
#   
#   # 准备残差
#   residual_c <- y_fit$residuals[(p1+1):length(obs_idx)]
#   residual_c <- residual_c - mean(residual_c)
#   
#   # 并行Bootstrap
#   Boot_m <- foreach(b = 1:B, .combine = "rbind") %dopar% {
#     boot_error <- sample(residual_c, Ti, replace = TRUE)
#     boot_y <- filter(boot_error, 
#                      filter = y_fit$coef[1:p1], 
#                      method = "recursive",
#                      init = rep(0, p1))
#     
#     for(t in (p1+1):Ti) {
#       x_term <- sum(y_fit$coef[-(1:p1)] * X_aug[t, ])
#       boot_y[t] <- boot_y[t] + x_term
#     }
#     
#     boot_fit <- arima(boot_y[obs_idx], order = c(p1,0,0), 
#                       xreg = X_aug_obs, include.mean = FALSE)
#     boot_pred <- predict(boot_fit, n.ahead = h, newxreg = X_future)$pred
#     boot_y[test_idx] - boot_pred
#   }
#   
#   # 计算置信区间
#   interval_m <- cbind(
#     predictions + apply(Boot_m, 2, quantile, probs = alpha/2, na.rm = TRUE),
#     predictions + apply(Boot_m, 2, quantile, probs = 1-alpha/2, na.rm = TRUE)
#   )
#   
#   list(predictions = predictions, interval = interval_m)
# }

Bootstrap_powered_ts <- function(X_m, y, tilde_y_m, p1, p2, obs_idx, test_idx,
                                 strong_idx, h, alpha=0.05, B=500){
  #X_m \in R^{T \times p}: input exogenous predictors
  #y \in R^T: the target time series
  #tilde_y_m \in R^{T \times q}: the synthetic surrogate
  #p1: the AR order for target ARX model
  #p2: the VAR order for synthetic surrogate
  #obs_idx: the observed index for target data
  #test_idx: the unobserved index for target data, the union of obs_idx and test_idx is equal to the whole index
  #strong_idx: the important index for predictors


  Ti <- dim(tilde_y_m)[1]
  invisible(capture.output(tilde_y_fit <- VARX(zt = tilde_y_m, p=p2, xt = X_m[,strong_idx], m=0, include.mean = FALSE)))
  ar_est <- tilde_y_fit$Phi
  hat_y <- matrix(0,nrow = Ti, dim(tilde_y_m)[2])
  for(t in (p2+1):Ti){
    hat_y[t,] <- (diag(ncol(tilde_y_m)) - ar_est)%*%t(tilde_y_m[t,,drop = FALSE])
  }

  y_obs <- y[obs_idx]
  X_aug <- cbind(X_m[,strong_idx], hat_y) # only use highly correlated index
  X_aug_obs <- X_aug[obs_idx,]

  y_fit <- arima(y_obs, order = c(p1,0,0), xreg = X_aug_obs, include.mean = FALSE)

  residual_c <- y_fit$residuals
  residual_c <- residual_c[(p1+1):length(obs_idx)]
  residual_c <- scale(residual_c, scale = FALSE) # demean

  X_future <- X_aug[test_idx,]
  predictions <- predict(y_fit, n.ahead = h, newxreg = X_future)$pred

  Boot_m <- matrix(nrow = B, ncol = h)
  for(b in 1:B){
    # bootstrap refit
    ## construct bootstrap dataset
    boot_error <- sample(residual_c, replace = TRUE, size = Ti)
    boot_y <- numeric(Ti)
    boot_y[1:p1] <- boot_error[1:p1]
    # regenerate ARX
    for (t in (p1 + 1):Ti) {
      # AR部分
      ar_term <- sum(y_fit$coef[1:p1] * boot_y[(t - 1):(t - p1)])

      # X部分
      x_term <- sum(y_fit$coef[(p1+1): length(y_fit$coef)] * X_aug[t, ])

      # 生成y_t
      boot_y[t] <- ar_term + x_term + boot_error[t]
    }

    # bootstrap estimation
    boot_y_fit <- arima(boot_y[obs_idx], order = c(p1,0,0), xreg = X_aug_obs, include.mean = FALSE)
    boot_predictions <- predict(boot_y_fit, n.ahead = h, newxreg = X_future)$pred

    # bootstrap error
    Boot_m[b,] <- boot_y[test_idx] - boot_predictions
   }


  interval_m <- matrix(nrow = h, ncol = 2)

  interval_m[,1] <- predictions + apply(Boot_m,2,quantile,alpha/2)
  interval_m[,2] <- predictions + apply(Boot_m,2,quantile,1-alpha/2)

  return(list(predictions, interval = interval_m))
}


# Bootstrap_powered_ts <- function(X_m, y, tilde_y_m, p1, p2, obs_idx, test_idx, 
#                                  strong_idx, h, alpha=0.05, B=500){
#   #X_m \in R^{T \times p}: input exogenous predictors
#   #y \in R^T: the target time series
#   #tilde_y_m \in R^{T \times q}: the synthetic surrogate 
#   #p1: the AR order for target ARX model
#   #p2: the VAR order for synthetic surrogate 
#   #obs_idx: the observed index for target data
#   #test_idx: the unobserved index for target data, the union of obs_idx and test_idx is equal to the whole index
#   #strong_idx: the important index for predictors
#   
#   
#   Ti <- dim(tilde_y_m)[1]
#   invisible(capture.output(tilde_y_fit <- VARX(zt = tilde_y_m, p=p2, xt = X_m[,strong_idx], m=0, include.mean = FALSE)))
#   
#   ar_est <- tilde_y_fit$Phi
#   hat_y <- matrix(0,nrow = Ti, dim(tilde_y_m)[2])
#   for(t in (p2+1):Ti){
#     hat_y[t,] <- (diag(ncol(tilde_y_m)) - ar_est)%*%t(tilde_y_m[t,,drop = FALSE])
#   }
#   
#   y_obs <- y[obs_idx]
#   X_aug <- cbind(X_m[,strong_idx], hat_y) # only use highly correlated index
#   X_aug_obs <- X_aug[obs_idx,]
#   
#   y_fit <- arima(y_obs, order = c(p1,0,0), xreg = X_aug_obs, include.mean = FALSE)
#   
#   residual_c <- y_fit$residuals 
#   residual_c <- residual_c[(p1+1):length(obs_idx)]
#   
#   X_future <- X_aug[test_idx,]
#   predictions <- predict(y_fit, n.ahead = h, newxreg = X_future)$pred
#   
#   Boot_m <- matrix(nrow = B, ncol = h)
#   for(b in 1:B){
#     res_residual <- sample(residual_c, h, replace = TRUE)
#     res_boot_pre_c <- c(y_obs[(length(y_obs)-p1+1):length(y_obs)],rep(0,h))
#     for(hh in 1:h){
#       res_ar <- 0
#       for(l in 1:p1){
#         res_ar <- res_ar + y_fit$coef[l]*res_boot_pre_c[hh+p1-l]
#       }
#       res_boot_pre_c[p1+hh] <- res_ar + sum(X_future[hh,]*y_fit$coef[(p1+1):length(y_fit$coef)]) + res_residual[hh]
#     }
#     Boot_m[b,] <- res_boot_pre_c[(p1+1):(h+p1)]
#   }
#   
#   interval_m <- matrix(nrow = h, ncol = 2)
#   
#   interval_m[,1] <- apply(Boot_m,2,quantile,alpha/2)
#   interval_m[,2] <- apply(Boot_m,2,quantile,1-alpha/2)
#   return(list(predictions, interval = interval_m))
# }

# RW_prediction <- function(y_obs, h, alpha = 0.05){
#   
#   sd <- sqrt(var(y_obs - lag(y_obs,1), na.rm = TRUE))
#   predictions <- rep(y_obs[length(y_obs)], h)
#   
#   interval_m <- matrix(nrow = h, ncol = 2)
#   var_h <- 0
#   for(i in 1:h){
#     interval_m[i,1] <- predictions[i]-sd*qnorm(alpha/2, lower.tail = FALSE)
#     interval_m[i,2] <- predictions[i]+sd*qnorm(alpha/2, lower.tail = FALSE)
#   }
#   return(list(predictions, interval = interval_m))
# }

Average_prediction <- function(y, H){
  pre_c <- rep(0, H)
  for(h in 1:H){
    pre_c[h] <- mean(y[length(y):(length(y)-h+1)])
  }
  return(pre_c)
} 


check_coverage <- function(pred_intervals, true_values) {
  coverage <- (true_values >= pred_intervals[1]) & (true_values <= pred_intervals[2])
  return(as.integer(coverage))
}

library(dplyr)
library(forecast)
library(MTS)
library(lubridate)
library(expm)
library(mvtnorm)
library(MASS)
