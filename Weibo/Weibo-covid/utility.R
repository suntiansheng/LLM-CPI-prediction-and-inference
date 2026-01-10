Prediction_powered_ts <- function(X_m, y, tilde_y_m, p1, p2, obs_idx, test_idx, 
                                  strong_idx, h){
  Ti <- dim(tilde_y_m)[1]



  invisible(capture.output(tilde_y_fit <- VARX(zt = tilde_y_m, p=p2, xt = X_m[,strong_idx], m=0, include.mean = FALSE)))

  

  
  ar_est <- tilde_y_fit$Phi
  hat_y <- matrix(0,nrow = Ti, dim(tilde_y_m)[2])
  for(t in (p2+1):Ti){
    hat_y[t,] <- (diag(ncol(tilde_y_m)) - ar_est)%*%t(tilde_y_m[t,,drop = FALSE])
  }
  
  y_obs <- y[obs_idx]
  X_aug <- cbind(X_m[,strong_idx], hat_y)
  X_aug_obs <- X_aug[obs_idx,]
  


  y_fit <- arima(y_obs, order = c(p1,0,0), xreg = X_aug_obs, include.mean = FALSE)
  

  


  X_future <- X_aug[test_idx,]

  
  predictions <- predict(y_fit, n.ahead = h, newxreg = X_future)$pred
  return(predictions)
}

Prediction_BJ <- function(X_m, y, tilde_y_m, p1, p2, obs_idx, test_idx, 
                          strong_idx, h, alpha=0.05){








  
  
  Ti <- dim(tilde_y_m)[1]


  invisible(capture.output(tilde_y_fit <- VARX(zt = tilde_y_m, p=p2, xt = X_m[,strong_idx], m=0, include.mean = FALSE)))

  

  
  ar_est <- tilde_y_fit$Phi
  hat_y <- matrix(0,nrow = Ti, dim(tilde_y_m)[2])
  for(t in (p2+1):Ti){
    hat_y[t,] <- (diag(ncol(tilde_y_m)) - ar_est)%*%t(tilde_y_m[t,,drop = FALSE])
  }
  
  y_obs <- y[obs_idx]
  X_aug <- cbind(X_m[,strong_idx], hat_y)
  X_aug_obs <- X_aug[obs_idx,]
  


  y_fit <- arima(y_obs, order = c(p1,0,0), xreg = X_aug_obs, include.mean = FALSE)
  
  var_e <- mean(y_fit$residuals^2)
  alpha_hat <- y_fit$coef[1:p1]
  I_q1 <- diag(p1-1)
  I_q1 <- cbind(I_q1, 0)
  A_hat <- rbind(alpha_hat, I_q1)
  X_future <- X_aug[test_idx,]

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


Bootstrap_powered_ts <- function(X_m, y, tilde_y_m, p1, p2, obs_idx, test_idx, 
                                 strong_idx, h, alpha=0.05, B=500){








  
  
  Ti <- dim(tilde_y_m)[1]
  invisible(capture.output(tilde_y_fit <- VARX(zt = tilde_y_m, p=p2, xt = X_m[,strong_idx], m=0, include.mean = FALSE)))
  
  ar_est <- tilde_y_fit$Phi
  hat_y <- matrix(0,nrow = Ti, dim(tilde_y_m)[2])
  for(t in (p2+1):Ti){
    hat_y[t,] <- (diag(ncol(tilde_y_m)) - ar_est)%*%t(tilde_y_m[t,,drop = FALSE])
  }
  
  y_obs <- y[obs_idx]
  X_aug <- cbind(X_m[,strong_idx], hat_y)
  X_aug_obs <- X_aug[obs_idx,]
  
  y_fit <- arima(y_obs, order = c(p1,0,0), xreg = X_aug_obs, include.mean = FALSE)
  
  residual_c <- y_fit$residuals 
  residual_c <- residual_c[(p1+1):length(obs_idx)]
  
  X_future <- X_aug[test_idx,]
  predictions <- predict(y_fit, n.ahead = h, newxreg = X_future)$pred
  
  Boot_m <- matrix(nrow = B, ncol = h)
  for(b in 1:B){
    res_residual <- sample(residual_c, h, replace = TRUE)
    res_boot_pre_c <- c(y_obs[(length(y_obs)-p1+1):length(y_obs)],rep(0,h))
    for(hh in 1:h){
      res_ar <- 0
      for(l in 1:p1){
        res_ar <- res_ar + y_fit$coef[l]*res_boot_pre_c[hh+p1-l]
      }
      res_boot_pre_c[p1+hh] <- res_ar + sum(X_future[hh,]*y_fit$coef[(p1+1):length(y_fit$coef)]) + res_residual[hh]
    }
    Boot_m[b,] <- res_boot_pre_c[(p1+1):(h+p1)]
  }
  
  interval_m <- matrix(nrow = h, ncol = 2)
  
  interval_m[,1] <- apply(Boot_m,2,quantile,alpha/2)
  interval_m[,2] <- apply(Boot_m,2,quantile,1-alpha/2)
  return(list(predictions, interval = interval_m))
}















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
