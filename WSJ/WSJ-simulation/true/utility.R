
suppressPackageStartupMessages({
  library(forecast)
  library(MTS)
  library(expm)
  library(dplyr)
})


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
  residual_c <- scale(residual_c, scale = FALSE)

  X_future <- X_aug[test_idx, , drop = FALSE]
  predictions <- predict(y_fit, n.ahead = h, newxreg = X_future)$pred

  Boot_m <- matrix(nrow = B, ncol = h)
  for (b in 1:B) {
    boot_error <- sample(residual_c, replace = TRUE, size = Ti)
    boot_y <- numeric(Ti)
    boot_y[1:p1] <- boot_error[1:p1]
    for (t in (p1 + 1):Ti) {
      ar_term <- sum(y_fit$coef[1:p1] * boot_y[(t - 1):(t - p1)])
      x_term <- sum(y_fit$coef[(p1 + 1):length(y_fit$coef)] * X_aug[t, ])
      boot_y[t] <- ar_term + x_term + boot_error[t]
    }
    boot_y_fit <- arima(boot_y[obs_idx], order = c(p1, 0, 0), xreg = X_aug_obs, include.mean = FALSE)
    boot_predictions <- predict(boot_y_fit, n.ahead = h, newxreg = X_future)$pred
    Boot_m[b, ] <- boot_y[test_idx] - boot_predictions
  }
  qs <- apply(Boot_m, 2, quantile, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
  interval_m <- cbind(predictions + qs[1, ], predictions + qs[2, ])
  list(predictions = predictions, interval = interval_m)
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

