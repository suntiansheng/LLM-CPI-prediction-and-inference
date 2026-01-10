
suppressPackageStartupMessages({
  library(forecast)
  library(MTS)
  library(dplyr)
})

if (!exists("%^%")) {
  `%^%` <- function(mat, n) {
    if (n < 0) {
      stop("Matrix power only supports non-negative integers.")
    }
    if (n == 0) {
      return(diag(nrow(mat)))
    }
    out <- mat
    if (n == 1) {
      return(out)
    }
    for (i in 2:n) {
      out <- out %*% mat
    }
    out
  }
}

build_hat_y <- function(tilde_y_m, X_m, strong_idx, p2) {
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
    acc <- Reduce(`+`, lapply(seq_len(p_use), function(l) {
      ylag <- as.numeric(tilde_y_m[t - l, ])
      Phi_list[[l]] %*% ylag
    }))
    hat_y[t, ] <- as.numeric(tilde_y_m[t, ]) - as.numeric(acc)
  }
  hat_y
}

build_arx_inputs <- function(X_m, y, tilde_y_m, p2, strong_idx, obs_idx) {
  hat_y <- build_hat_y(tilde_y_m, X_m, strong_idx, p2)
  X_aug <- cbind(X_m[, strong_idx, drop = FALSE], hat_y)
  list(
    hat_y = hat_y,
    X_aug = X_aug,
    X_aug_obs = X_aug[obs_idx, , drop = FALSE],
    y_obs = y[obs_idx]
  )
}

recursive_arx_forecast <- function(y, X, coef_ar, coef_x, p, start_idx, h) {
  if (start_idx + h > nrow(X)) {
    stop("Not enough future surrogate rows in X for recursive forecast.")
  }
  y_ext <- y[1:start_idx]
  yhat <- numeric(h)
  for (k in seq_len(h)) {
    lags <- rev(tail(y_ext, p))
    ar_term <- sum(coef_ar * lags)
    x_term <- sum(coef_x * X[start_idx + k, ])
    yhat[k] <- ar_term + x_term
    y_ext <- c(y_ext, yhat[k])
  }
  yhat
}


LLM_TS.Predict <- function(X_m, y, tilde_y_m, p1, p2, obs_idx, test_idx,
                                  strong_idx, h) {
  arx_inputs <- build_arx_inputs(X_m, y, tilde_y_m, p2, strong_idx, obs_idx)
  y_fit <- arima(arx_inputs$y_obs, order = c(p1, 0, 0), xreg = arx_inputs$X_aug_obs, include.mean = FALSE)
  X_future <- arx_inputs$X_aug[test_idx, , drop = FALSE]
  predict(y_fit, n.ahead = h, newxreg = X_future)$pred
}

LLM_TS.BJ <- function(X_m, y, tilde_y_m, p1, p2, obs_idx, test_idx,
                      strong_idx, h, alpha = 0.05) {
  arx_inputs <- build_arx_inputs(X_m, y, tilde_y_m, p2, strong_idx, obs_idx)
  y_fit <- arima(arx_inputs$y_obs, order = c(p1,0,0), xreg = arx_inputs$X_aug_obs, include.mean = FALSE)
  coef <- y_fit$coef
  coef_ar <- coef[1:p1]
  coef_x <- coef[(p1+1):length(coef)]


  n_obs <- length(obs_idx)
  n_cal <- max(5, floor(0.15 * n_obs))
  n_fit <- n_obs - n_cal
  if (n_fit <= p1 + 1) {
    stop("Not enough data for BJ calibration split.")
  }
  cal_start <- n_fit + 1
  var_h_vec <- rep(NA_real_, h)
  for (i in seq_len(h)) {
    max_origin <- n_obs - i
    if (max_origin < cal_start) next
    origins <- obs_idx[cal_start:max_origin]
    resids <- vapply(origins, function(origin) {
      pred_i <- recursive_arx_forecast(y, arx_inputs$X_aug, coef_ar, coef_x, p1, origin, i)
      y[origin + i] - pred_i[i]
    }, numeric(1))
    var_h_vec[i] <- mean(resids^2, na.rm = TRUE)
  }



  start_idx <- max(obs_idx)

  predictions <- recursive_arx_forecast(
    y = y,
    X = arx_inputs$X_aug,
    coef_ar = coef_ar,
    coef_x = coef_x,
    p = p1,
    start_idx = start_idx,
    h = h
  )


  interval_m <- matrix(nrow = h, ncol = 2)
  fallback_var <- mean(y_fit$residuals^2, na.rm = TRUE)
  for (i in seq_len(h)) {
    var_h <- if (is.na(var_h_vec[i])) fallback_var else var_h_vec[i]
    n_resid_i <- max(0, (n_obs - i) - cal_start + 1)
    df_i <- if (n_resid_i > 2) n_resid_i - 1 else max(2, length(y_fit$residuals) - 1)
    t_crit <- qt(1 - alpha / 2, df = df_i)
    sd_h <- sqrt(var_h)
    interval_m[i, 1] <- predictions[i] - t_crit * sd_h
    interval_m[i, 2] <- predictions[i] + t_crit * sd_h
  }

  list(predictions = predictions, interval = interval_m)
}






















































LLM_TS.BOOT <- function(X_m, y, tilde_y_m, p1, p2, obs_idx, test_idx,
                                 strong_idx, h, alpha = 0.05, B = 500) {
  Ti <- nrow(tilde_y_m)
  arx_inputs <- build_arx_inputs(X_m, y, tilde_y_m, p2, strong_idx, obs_idx)
  y_fit <- arima(arx_inputs$y_obs, order = c(p1, 0, 0), xreg = arx_inputs$X_aug_obs, include.mean = FALSE)

  residual_c <- y_fit$residuals
  residual_c <- residual_c[(p1 + 1):length(obs_idx)]

  X_future <- arx_inputs$X_aug[test_idx, , drop = FALSE]
  predictions <- predict(y_fit, n.ahead = h, newxreg = X_future)$pred

  Boot_m <- matrix(nrow = B, ncol = h)
  boot_y_init <- y[1:p1]
  for (b in 1:B) {
    boot_error <- sample(residual_c, replace = TRUE, size = Ti)
    boot_y <- numeric(Ti)
    boot_y[1:p1] <- boot_y_init
    for (t in (p1 + 1):Ti) {
      ar_term <- sum(y_fit$coef[1:p1] * boot_y[(t - 1):(t - p1)])
      x_term <- sum(y_fit$coef[(p1 + 1):length(y_fit$coef)] * arx_inputs$X_aug[t, ])
      boot_y[t] <- ar_term + x_term + boot_error[t]
    }
    boot_y_fit <- arima(boot_y[obs_idx], order = c(p1, 0, 0), xreg = arx_inputs$X_aug_obs, include.mean = FALSE)
    boot_predictions <- predict(boot_y_fit, n.ahead = h, newxreg = X_future)$pred
    Boot_m[b, ] <- boot_predictions - boot_y[test_idx]
  }
  qs <- apply(Boot_m, 2, quantile, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
  interval_m <- rbind(predictions - qs[2, ], predictions - qs[1, ])
  list(predictions = predictions, interval = t(interval_m))
}




select_topics <- function(y, X, train_idx, use_cor = TRUE,
                          max_topics = Inf, penalty_w = 0.1) {

  ar_base <- auto.arima(y[train_idx], max.q = 0, D = 0, seasonal = FALSE,
                        allowmean = FALSE, allowdrift = FALSE)
  p_ar <- ar_base$arma[1]
  ar_fit <- arima(y[train_idx], order = c(p_ar, 0, 0), include.mean = FALSE)
  resid_base <- residuals(ar_fit)

  scores <- sapply(colnames(X), function(col) {
    z <- X[train_idx, col]
    if (use_cor) cor(resid_base, z, use = "complete.obs") else cov(resid_base, z, use = "complete.obs")
  })
  topic_order <- names(sort(abs(scores), decreasing = TRUE))

  selected <- character()
  best_aic <- Inf
  path <- list()

  for (cand in topic_order) {
    if (length(selected) >= max_topics) break
    new_set <- c(selected, cand)
    xi_new <- X[train_idx, new_set, drop = FALSE]
    fit_new <- tryCatch(
      arima(y[train_idx], order = c(p_ar, 0, 0), xreg = xi_new, include.mean = FALSE),
      error = function(e) NULL
    )
    if (is.null(fit_new)) next
    m <- length(new_set)
    t1 <- length(train_idx)
    pen <- penalty_w * ((m + 1) * (m + 2)) / (t1 - m - 2)
    cand_aic <- fit_new$aic + pen
    path[[length(path) + 1]] <- data.frame(
      step = length(new_set),
      topic_added = cand,
      combo = paste(new_set, collapse = ","),
      aic_pen = cand_aic,
      raw_aic = fit_new$aic,
      score = scores[cand],
      stringsAsFactors = FALSE
    )
    if (cand_aic < best_aic) {
      selected <- new_set
      best_aic <- cand_aic
    } else {
      break
    }
  }

  list(selected = selected,
       ar_order = p_ar,
       path = if (length(path)) do.call(rbind, path) else data.frame())
}
