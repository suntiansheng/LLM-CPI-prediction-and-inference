#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(forecast)
  library(MTS)
  library(expm)
  library(dplyr)
})

safe_auto_arima <- function(...) {
  fit <- try(suppressWarnings(auto.arima(...)), silent = TRUE)
  if (inherits(fit, "try-error")) {
    return(list(arma = c(1, 0, 0)))
  }
  fit
}

LLM_TS.Predict <- function(X_m, y, tilde_y_m, p1, p2, obs_idx, test_idx, 
                                  strong_idx, h) {
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

LLM_TS.BJ <- function(X_m, y, tilde_y_m, p1, p2, obs_idx, test_idx, 
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

LLM_TS.BOOT <- function(X_m, y, tilde_y_m, p1, p2, obs_idx, test_idx,
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
    boot_error <- sample(residual_c, replace = TRUE, size = Ti)
    boot_y <- numeric(Ti)
    for (t in (p1 + 1):Ti) {
      ar_term <- sum(y_fit$coef[1:p1] * boot_y[(t - 1):(t - p1)])
      x_term <- sum(y_fit$coef[(p1+1): length(y_fit$coef)] * X_aug[t, ])
      boot_y[t] <- ar_term + x_term + boot_error[t]
    }

    boot_y_fit <- arima(boot_y[obs_idx], order = c(p1,0,0), xreg = X_aug_obs, include.mean = FALSE)
    boot_predictions <- predict(boot_y_fit, n.ahead = h, newxreg = X_future)$pred

    Boot_m[b,] <- boot_predictions - boot_y[test_idx]
   }

  interval_m <- matrix(nrow = h, ncol = 2)
  interval_m[,1] <- predictions + apply(Boot_m,2,quantile,alpha/2)
  interval_m[,2] <- predictions + apply(Boot_m,2,quantile,1-alpha/2)

  return(list(predictions, interval = interval_m))
}

LLM_TS.Predict <- function(X_m, y, tilde_y_m, p1, p2, obs_idx, test_idx,
                                  strong_idx, h) {
  Ti <- nrow(tilde_y_m)
  invisible(capture.output(
    tilde_y_fit <- VARX(zt = tilde_y_m, p = p2, xt = X_m[, strong_idx, drop = FALSE],
                        m = 0, include.mean = FALSE)
  ))
  Phi_arr <- tilde_y_fit$Phi
  Phi_list <- if (length(dim(Phi_arr)) == 3) {
    lapply(seq_len(dim(Phi_arr)[3]), function(j) as.matrix(Phi_arr[, , j]))
  } else {
    list(as.matrix(Phi_arr))
  }
  p_use <- length(Phi_list)
  hat_y <- matrix(0, nrow = Ti, ncol = ncol(tilde_y_m))
  for (t in seq_len(Ti)) {
    if (t <= p_use) next
    acc <- Reduce(`+`, lapply(seq_len(p_use), function(l) Phi_list[[l]] %*% tilde_y_m[t - l, ]))
    hat_y[t, ] <- tilde_y_m[t, ] - acc
  }

  y_obs <- y[obs_idx]
  X_aug <- cbind(X_m[, strong_idx, drop = FALSE], hat_y)
  X_aug_obs <- X_aug[obs_idx, , drop = FALSE]

  y_fit <- arima(y_obs, order = c(p1, 0, 0), xreg = X_aug_obs, include.mean = FALSE)
  X_future <- X_aug[test_idx, , drop = FALSE]
  predict(y_fit, n.ahead = h, newxreg = X_future)$pred
}

LLM_TS.BJ <- function(X_m, y, tilde_y_m, p1, p2, obs_idx, test_idx,
                          strong_idx, h, alpha = 0.05) {
  Ti <- nrow(tilde_y_m)
  invisible(capture.output(
    tilde_y_fit <- VARX(zt = tilde_y_m, p = p2, xt = X_m[, strong_idx, drop = FALSE],
                        m = 0, include.mean = FALSE)
  ))

  Phi_arr <- tilde_y_fit$Phi
  Phi_list <- if (length(dim(Phi_arr)) == 3) {
    lapply(seq_len(dim(Phi_arr)[3]), function(j) as.matrix(Phi_arr[, , j]))
  } else {
    list(as.matrix(Phi_arr))
  }
  p_use <- length(Phi_list)
  hat_y <- matrix(0, nrow = Ti, ncol = ncol(tilde_y_m))
  for (t in seq_len(Ti)) {
    if (t <= p_use) next
    acc <- Reduce(`+`, lapply(seq_len(p_use), function(l) Phi_list[[l]] %*% tilde_y_m[t - l, ]))
    hat_y[t, ] <- tilde_y_m[t, ] - acc
  }

  y_obs <- y[obs_idx]
  X_aug <- cbind(X_m[, strong_idx, drop = FALSE], hat_y)
  X_aug_obs <- X_aug[obs_idx, , drop = FALSE]

  y_fit <- arima(y_obs, order = c(p1, 0, 0), xreg = X_aug_obs, include.mean = FALSE)

  var_e <- mean(y_fit$residuals^2)
  alpha_hat <- y_fit$coef[1:p1]
  I_q1 <- diag(p1 - 1)
  I_q1 <- cbind(I_q1, 0)
  A_hat <- rbind(alpha_hat, I_q1)
  X_future <- X_aug[test_idx, , drop = FALSE]
  predictions <- predict(y_fit, n.ahead = h, newxreg = X_future)$pred

  interval_m <- matrix(nrow = h, ncol = 2)
  var_h <- 0
  for (i in 1:h) {
    res_var <- var_e * (A_hat %^% (i - 1))[1, 1]
    var_h <- res_var + var_h
    interval_m[i, 1] <- predictions[i] - sqrt(var_h) * qnorm(alpha / 2, lower.tail = FALSE)
    interval_m[i, 2] <- predictions[i] + sqrt(var_h) * qnorm(alpha / 2, lower.tail = FALSE)
  }
  list(predictions = predictions, interval = interval_m)
}

LLM_TS.BOOT <- function(X_m, y, tilde_y_m, p1, p2, obs_idx, test_idx,
                                 strong_idx, h, alpha = 0.05, B = 500) {
  Ti <- nrow(tilde_y_m)
  invisible(capture.output(
    tilde_y_fit <- VARX(zt = tilde_y_m, p = p2, xt = X_m[, strong_idx, drop = FALSE],
                        m = 0, include.mean = FALSE)
  ))
  Phi_arr <- tilde_y_fit$Phi
  Phi_list <- if (length(dim(Phi_arr)) == 3) {
    lapply(seq_len(dim(Phi_arr)[3]), function(j) as.matrix(Phi_arr[, , j]))
  } else {
    list(as.matrix(Phi_arr))
  }
  p_use <- length(Phi_list)
  hat_y <- matrix(0, nrow = Ti, ncol = ncol(tilde_y_m))
  for (t in seq_len(Ti)) {
    if (t <= p_use) next
    acc <- Reduce(`+`, lapply(seq_len(p_use), function(l) Phi_list[[l]] %*% tilde_y_m[t - l, ]))
    hat_y[t, ] <- tilde_y_m[t, ] - acc
  }

  y_obs <- y[obs_idx]
  X_aug <- cbind(X_m[, strong_idx, drop = FALSE], hat_y)
  X_aug_obs <- X_aug[obs_idx, , drop = FALSE]

  y_fit <- arima(y_obs, order = c(p1, 0, 0), xreg = X_aug_obs, include.mean = FALSE)

  residual_c <- y_fit$residuals
  residual_c <- residual_c[(p1 + 1):length(obs_idx)]

  X_future <- X_aug[test_idx, , drop = FALSE]
  predictions <- predict(y_fit, n.ahead = h, newxreg = X_future)$pred

  Boot_m <- matrix(nrow = B, ncol = h)
  boot_y_init <- y[1:p1]
  for (b in 1:B) {
    boot_error <- sample(residual_c, replace = TRUE, size = Ti)
    boot_y <- numeric(Ti)
    boot_y[1:p1] <- boot_y_init
    for (t in (p1 + 1):Ti) {
      ar_term <- sum(y_fit$coef[1:p1] * boot_y[(t - 1):(t - p1)])
      x_term <- sum(y_fit$coef[(p1 + 1):length(y_fit$coef)] * X_aug[t, ])
      boot_y[t] <- ar_term + x_term + boot_error[t]
    }
    boot_y_fit <- arima(boot_y[obs_idx], order = c(p1, 0, 0), xreg = X_aug_obs, include.mean = FALSE)
    boot_predictions <- predict(boot_y_fit, n.ahead = h, newxreg = X_future)$pred
    Boot_m[b, ] <- boot_predictions - boot_y[test_idx]
  }
  qs <- apply(Boot_m, 2, quantile, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
  interval_m <- rbind(predictions - qs[2, ], predictions - qs[1, ])
  list(predictions = predictions, interval = t(interval_m))
}

Average_prediction <- function(y, H) {
  sapply(seq_len(H), function(h) mean(tail(y, h)))
}

check_coverage <- function(interval_row, true_val) {
  as.integer(true_val >= interval_row[1] & true_val <= interval_row[2])
}

run_prediction_eval <- function(y, X_m, unem, tilde_y_m, strong_idx, h_vec,
                                has_unem = TRUE, p2 = 1) {
  y_vec <- as.numeric(y)
  unem_vec <- if (has_unem) as.numeric(unem) else NULL
  Res_m <- matrix(NA_real_, nrow = length(h_vec), ncol = 16)

  for (k in seq_along(h_vec)) {
    h <- h_vec[k]
    res_c <- NULL
    test_idx <- (length(y_vec) - h + 1):length(y_vec)
    obs_idx <- setdiff(seq_along(y_vec), test_idx)

    ar_base <- tryCatch(auto.arima(y_vec[obs_idx], max.q = 0, D = 0, seasonal = FALSE,
                                   allowmean = FALSE, allowdrift = FALSE),
                        error = function(e) list(arma = c(1, 0, 0)))
    ar_fit <- arima(y_vec[obs_idx], order = c(ar_base$arma[1], 0, 0), include.mean = FALSE)
    ar_prediction <- predict(ar_fit, n.ahead = h)$pred
    res_c <- c(res_c, sqrt(mean((ar_prediction - y_vec[test_idx])^2)))
    res_c <- c(res_c, mean(sign(ar_prediction) != sign(y_vec[test_idx])))

    mean_prediction <- rep(y_vec[max(obs_idx)], length(test_idx))
    res_c <- c(res_c, sqrt(mean((mean_prediction - y_vec[test_idx])^2)))
    res_c <- c(res_c, mean(sign(mean_prediction) != sign(y_vec[test_idx])))

    ave_prediction <- Average_prediction(y_vec[obs_idx], length(test_idx))
    res_c <- c(res_c, sqrt(mean((ave_prediction - y_vec[test_idx])^2)))
    res_c <- c(res_c, mean(sign(ave_prediction) != sign(y_vec[test_idx])))

    if (has_unem) {
      ar_unem_base <- tryCatch(auto.arima(y_vec[obs_idx], xreg = unem_vec[obs_idx], max.q = 0,
                                          D = 0, seasonal = FALSE, allowmean = FALSE,
                                          allowdrift = FALSE),
                               error = function(e) list(arma = c(1, 0, 0)))
      ar_unem_fit <- arima(y_vec[obs_idx], xreg = unem_vec[obs_idx],
                           order = c(ar_unem_base$arma[1], 0, 0), include.mean = FALSE)
      ar_unem_prediction <- predict(ar_unem_fit, newxreg = unem_vec[test_idx], n.ahead = h)$pred
      res_c <- c(res_c, sqrt(mean((ar_unem_prediction - y_vec[test_idx])^2)))
      res_c <- c(res_c, mean(sign(ar_unem_prediction) != sign(y_vec[test_idx])))
    } else {
      res_c <- c(res_c, NA_real_, NA_real_)
    }

    arx_fit <- arima(y_vec[obs_idx], order = c(ar_base$arma[1], 0, 0),
                     xreg = X_m[obs_idx, strong_idx, drop = FALSE], include.mean = FALSE)
    arx_predictions <- predict(arx_fit, n.ahead = h,
                               newxreg = X_m[test_idx, strong_idx, drop = FALSE])$pred
    res_c <- c(res_c, sqrt(mean((arx_predictions - y_vec[test_idx])^2)))
    res_c <- c(res_c, mean(sign(arx_predictions) != sign(y_vec[test_idx])))

    if (has_unem) {
      arx_unem_fit <- arima(y_vec[obs_idx], order = c(max(1, ar_base$arma[1]), 0, 0),
                            xreg = cbind(X_m[obs_idx, strong_idx, drop = FALSE],
                                         unem_vec[obs_idx]),
                            include.mean = FALSE)
      arx_unem_predictions <- predict(arx_unem_fit, n.ahead = h,
                                      newxreg = cbind(X_m[test_idx, strong_idx, drop = FALSE],
                                                      unem_vec[test_idx]))$pred
      res_c <- c(res_c, sqrt(mean((arx_unem_predictions - y_vec[test_idx])^2)))
      res_c <- c(res_c, mean(sign(arx_unem_predictions) != sign(y_vec[test_idx])))
    } else {
      res_c <- c(res_c, NA_real_, NA_real_)
    }

    p1 <- ar_base$arma[1]
    powered_prediction <- LLM_TS.Predict(X_m, y_vec, tilde_y_m, p1, p2,
                                         obs_idx, test_idx, strong_idx, h)
    res_c <- c(res_c, sqrt(mean((powered_prediction - y_vec[test_idx])^2)))
    res_c <- c(res_c, mean(sign(powered_prediction) != sign(y_vec[test_idx])))

    if (has_unem) {
      powered_unem_prediction <- LLM_TS.Predict(cbind(X_m, unem_vec), y_vec, tilde_y_m,
                                                p1, p2, obs_idx, test_idx,
                                                c(strong_idx, ncol(X_m) + 1), h)
      res_c <- c(res_c, sqrt(mean((powered_unem_prediction - y_vec[test_idx])^2)))
      res_c <- c(res_c, mean(sign(powered_unem_prediction) != sign(y_vec[test_idx])))
    } else {
      res_c <- c(res_c, NA_real_, NA_real_)
    }

    Res_m[k, ] <- res_c
  }

  list(metrics = Res_m)
}

run_inference_eval <- function(y, X_m, unem, tilde_y_m, strong_idx, h_vec,
                               has_unem = TRUE, p2 = 1, alpha = 0.05, boot_B = 2000) {
  y_vec <- as.numeric(y)
  unem_vec <- if (has_unem) as.numeric(unem) else NULL
  coverage_m <- matrix(NA_real_, nrow = length(h_vec), ncol = 8)
  len_m <- matrix(NA_real_, nrow = length(h_vec), ncol = 8)

  for (k in seq_along(h_vec)) {
    h <- h_vec[k]
    test_idx <- (length(y_vec) - h + 1):length(y_vec)
    obs_idx <- setdiff(seq_along(y_vec), test_idx)

    ar_base <- tryCatch(auto.arima(y_vec[obs_idx], max.q = 0, D = 0, seasonal = FALSE,
                                   allowmean = FALSE, allowdrift = FALSE),
                        error = function(e) list(arma = c(1, 0, 0)))
    ar_fit <- Arima(y_vec[obs_idx], order = c(ar_base$arma[1], 0, 0), include.mean = FALSE)
    ar_predictions_fit <- forecast(ar_fit, h = length(test_idx))
    ar_interval <- cbind(ar_predictions_fit$lower[, 2], ar_predictions_fit$upper[, 2])
    len_m[k, 1] <- mean(as.numeric(ar_interval[, 2] - ar_interval[, 1]))
    coverage_m[k, 1] <- mean(mapply(check_coverage, split(ar_interval, row(ar_interval)),
                                    y_vec[test_idx]))

    if (has_unem) {
      ar_unem_fit <- Arima(y_vec[obs_idx], order = c(ar_base$arma[1], 0, 0),
                           xreg = unem_vec[obs_idx], include.mean = FALSE)
      ar_unem_predictions_fit <- forecast(ar_unem_fit, h = length(test_idx),
                                           xreg = unem_vec[test_idx])
      ar_unem_interval <- cbind(ar_unem_predictions_fit$lower[, 2],
                                ar_unem_predictions_fit$upper[, 2])
      len_m[k, 2] <- mean(as.numeric(ar_unem_interval[, 2] - ar_unem_interval[, 1]))
      coverage_m[k, 2] <- mean(mapply(check_coverage, split(ar_unem_interval, row(ar_unem_interval)),
                                      y_vec[test_idx]))
    }

    arx_lda_fit <- Arima(y_vec[obs_idx], order = c(ar_base$arma[1], 0, 0),
                         xreg = X_m[obs_idx, strong_idx, drop = FALSE], include.mean = FALSE)
    arx_lda_predictions_fit <- forecast(arx_lda_fit, h = length(test_idx),
                                        xreg = X_m[test_idx, strong_idx, drop = FALSE])
    arx_lda_interval <- cbind(arx_lda_predictions_fit$lower[, 2],
                              arx_lda_predictions_fit$upper[, 2])
    len_m[k, 3] <- mean(arx_lda_interval[, 2] - arx_lda_interval[, 1])
    coverage_m[k, 3] <- mean(mapply(check_coverage, split(arx_lda_interval, row(arx_lda_interval)),
                                    y_vec[test_idx]))

    if (has_unem) {
      arx_lda_un_fit <- Arima(y_vec[obs_idx], order = c(ar_base$arma[1], 0, 0),
                              xreg = cbind(X_m[obs_idx, strong_idx, drop = FALSE],
                                           unem_vec[obs_idx]),
                              include.mean = FALSE)
      arx_lda_un_predictions_fit <- forecast(arx_lda_un_fit, h = length(test_idx),
                                             xreg = cbind(X_m[test_idx, strong_idx, drop = FALSE],
                                                          unem_vec[test_idx]))
      arx_lda_un_interval <- cbind(arx_lda_un_predictions_fit$lower[, 2],
                                   arx_lda_un_predictions_fit$upper[, 2])
      len_m[k, 4] <- mean(arx_lda_un_interval[, 2] - arx_lda_un_interval[, 1])
      coverage_m[k, 4] <- mean(mapply(check_coverage, split(arx_lda_un_interval, row(arx_lda_un_interval)),
                                      y_vec[test_idx]))
    }

    p1 <- ar_base$arma[1]
    powered_lda_interval <- LLM_TS.BJ(X_m, y_vec, tilde_y_m, p1, p2,
                                      obs_idx, test_idx, strong_idx, h, alpha = alpha)$interval
    len_m[k, 5] <- mean(powered_lda_interval[, 2] - powered_lda_interval[, 1])
    coverage_m[k, 5] <- mean(mapply(check_coverage, split(powered_lda_interval, row(powered_lda_interval)),
                                    y_vec[test_idx]))

    boot_lda_interval <- LLM_TS.BOOT(X_m, y_vec, tilde_y_m, p1, p2,
                                    obs_idx, test_idx, strong_idx, h,
                                    alpha = alpha, B = boot_B)$interval
    len_m[k, 6] <- mean(boot_lda_interval[, 2] - boot_lda_interval[, 1])
    coverage_m[k, 6] <- mean(mapply(check_coverage, split(boot_lda_interval, row(boot_lda_interval)),
                                    y_vec[test_idx]))

    if (has_unem) {
      powered_interval <- LLM_TS.BJ(cbind(X_m, unem_vec), y_vec, tilde_y_m, p1, p2,
                                    obs_idx, test_idx, c(strong_idx, ncol(X_m) + 1), h,
                                    alpha = alpha)$interval
      len_m[k, 7] <- mean(powered_interval[, 2] - powered_interval[, 1])
      coverage_m[k, 7] <- mean(mapply(check_coverage, split(powered_interval, row(powered_interval)),
                                      y_vec[test_idx]))

      boot_interval <- LLM_TS.BOOT(cbind(X_m, unem_vec), y_vec, tilde_y_m, p1, p2,
                                   obs_idx, test_idx, c(strong_idx, ncol(X_m) + 1), h,
                                   alpha = alpha, B = boot_B)$interval
      len_m[k, 8] <- mean(boot_interval[, 2] - boot_interval[, 1])
      coverage_m[k, 8] <- mean(mapply(check_coverage, split(boot_interval, row(boot_interval)),
                                      y_vec[test_idx]))
    }
  }

  list(coverage = coverage_m, length = len_m)
}

select_topics <- function(y, X_m, EOD, add_one = TRUE) {
  ar_fit <- auto.arima(y[1:EOD], max.q = 0, D = 0, seasonal = FALSE,
                       allowmean = FALSE, allowdrift = FALSE)
  ar_fit <- arima(y[1:EOD], order = c(ar_fit$arma[1], 0, 0), include.mean = FALSE)
  aic_c <- NULL
  for (i in 1:20) {
    strong_idx <- order(abs(cov(ar_fit$residuals, X_m[1:EOD, ])), decreasing = TRUE)[1:i]
    arx_fit <- arima(y[1:EOD, ], order = c(ar_fit$arma[1], 0, 0),
                     xreg = X_m[1:EOD, strong_idx], include.mean = FALSE)
    aic_c[i] <- EOD * log(mean((arx_fit$residuals)^2) + 1) +
      (2 * (length(strong_idx) + 1) * (length(strong_idx) + 2)) /
      (EOD - length(strong_idx) - 2)
  }
  select_num <- which(diff(aic_c) < 0)[1]
  if (add_one) {
    select_num <- select_num + 1
  }
  strong_idx <- order(abs(cov(ar_fit$residuals, X_m[1:EOD, ])), decreasing = TRUE)[1:select_num]
  strong_idx
}
