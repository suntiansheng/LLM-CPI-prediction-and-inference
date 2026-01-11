












spi_transport_scores <- function(real_scores, synth_scores) {
  real_scores  <- as.numeric(real_scores)
  synth_scores <- as.numeric(synth_scores)
  
  n_real  <- length(real_scores)
  n_synth <- length(synth_scores)
  
  if (n_real < 1L || n_synth < 1L) {
    stop("Need at least one real and one synthetic score.")
  }
  

  F_real <- ecdf(real_scores)
  u <- F_real(real_scores)
  

  eps <- 1e-8
  u[u <= 0] <- eps
  u[u >= 1] <- 1 - eps
  


  z <- as.numeric(quantile(synth_scores, probs = u, type = 1, names = FALSE))
  return(z)
}


spi_threshold <- function(real_scores, synth_scores, alpha = 0.1) {
  real_scores  <- as.numeric(real_scores)
  synth_scores <- as.numeric(synth_scores)
  
  if (alpha <= 0 || alpha >= 1) {
    stop("alpha must be in (0,1).")
  }
  

  z_real <- spi_transport_scores(real_scores, synth_scores)
  

  q_synth <- as.numeric(quantile(synth_scores, probs = 1 - alpha,
                                 type = 1, names = FALSE))
  

  idx <- which(z_real <= q_synth + 1e-12)
  
  if (length(idx) == 0L) {


    q_spi <- min(real_scores)
  } else {

    q_spi <- max(real_scores[idx])
  }
  
  list(
    q_spi   = q_spi,
    q_synth = q_synth,
    z_real  = z_real
  )
}












spi_interval_one_step <- function(y_cal,
                                  yhat_cal,
                                  synth_scores,
                                  yhat_test,
                                  alpha = 0.1) {
  y_cal    <- as.numeric(y_cal)
  yhat_cal <- as.numeric(yhat_cal)
  
  if (length(y_cal) != length(yhat_cal)) {
    stop("y_cal and yhat_cal must have the same length.")
  }
  

  real_scores <- abs(y_cal - yhat_cal)
  
  thr <- spi_threshold(real_scores = real_scores,
                       synth_scores = synth_scores,
                       alpha = alpha)
  
  q_spi <- thr$q_spi
  

  lower <- yhat_test - q_spi
  upper <- yhat_test + q_spi
  
  list(
    lower      = lower,
    upper      = upper,
    radius_spi = q_spi,
    q_synth    = thr$q_synth,
    real_scores = real_scores,
    z_real     = thr$z_real
  )
}
































spi_multi_horizon <- function(y,
                              X_m,
                              tilde_y_m,
                              strong_idx,
                              h_vec,
                              unem = NULL,
                              has_unem = FALSE,
                              cal_frac = 0.2,
                              alpha = 0.05,
                              p2 = 1) {
  y <- as.numeric(y)
  X_m <- as.matrix(X_m)
  tilde_y_m <- as.matrix(tilde_y_m)
  n <- length(y)
  if (nrow(X_m) != n || nrow(tilde_y_m) != n) {
    stop("X_m, tilde_y_m, and y must have the same number of rows.")
  }
  if (has_unem && is.null(unem)) {
    stop("Unemployment series required when has_unem = TRUE.")
  }
  if (has_unem) unem <- as.numeric(unem)


  x_varx <- X_m[, strong_idx, drop = FALSE]
  if (has_unem) {
    x_varx <- cbind(x_varx, unem)
  }
  invisible(capture.output(
    tilde_y_fit <- VARX(zt = tilde_y_m, p = p2, xt = x_varx,
                        m = 0, include.mean = FALSE)
  ))
  Phi_arr <- tilde_y_fit$Phi
  Phi_list <- if (length(dim(Phi_arr)) == 3) {
    lapply(seq_len(dim(Phi_arr)[3]), function(j) as.matrix(Phi_arr[, , j]))
  } else {
    list(as.matrix(Phi_arr))
  }
  p_use <- length(Phi_list)
  hat_y <- matrix(0, nrow = n, ncol = ncol(tilde_y_m))
  for (t in seq_len(n)) {
    if (t <= p_use) next
    acc <- Reduce(`+`, lapply(seq_len(p_use), function(l) Phi_list[[l]] %*% tilde_y_m[t - l, ]))
    hat_y[t, ] <- tilde_y_m[t, ] - acc
  }

  X_aug <- cbind(X_m[, strong_idx, drop = FALSE], hat_y)
  X_aug_unem <- if (has_unem) cbind(X_m[, strong_idx, drop = FALSE], unem, hat_y) else NULL

  cov_no <- len_no <- rep(NA_real_, length(h_vec))
  cov_un <- len_un <- rep(NA_real_, length(h_vec))

  for (k in seq_along(h_vec)) {
    h <- h_vec[k]
    test_idx <- (n - h + 1):n
    obs_idx <- setdiff(seq_len(n), test_idx)

    n_obs <- length(obs_idx)
    n_cal <- max(ceiling(cal_frac * n_obs), h)
    n_train <- n_obs - n_cal
    if (n_train <= 1) next

    train_idx <- obs_idx[seq_len(n_train)]
    cal_idx <- obs_idx[(n_train + 1):n_obs]


    ar_base <- tryCatch(
      auto.arima(y[train_idx], max.q = 0, D = 0, seasonal = FALSE,
                 allowmean = FALSE, allowdrift = FALSE),
      error = function(e) list(arma = c(1, 0, 0))
    )
    p1 <- ar_base$arma[1]
    fit_no <- tryCatch(
      arima(y[train_idx], order = c(p1, 0, 0),
            xreg = X_aug[train_idx, , drop = FALSE],
            include.mean = FALSE),
      error = function(e) NULL
    )
    if (!is.null(fit_no)) {
      yhat_cal <- predict(fit_no, n.ahead = length(cal_idx),
                          newxreg = X_aug[cal_idx, , drop = FALSE])$pred
      real_scores <- abs(y[cal_idx] - yhat_cal)
      synth_scores <- as.numeric(abs(hat_y[cal_idx, , drop = FALSE]))
      thr_no <- spi_threshold(real_scores = real_scores,
                              synth_scores = synth_scores,
                              alpha = alpha)
      yhat_test <- predict(fit_no, n.ahead = h,
                           newxreg = X_aug[test_idx, , drop = FALSE])$pred
      spi_interval <- cbind(yhat_test - thr_no$q_spi, yhat_test + thr_no$q_spi)
      len_no[k] <- mean(spi_interval[, 2] - spi_interval[, 1])
      cov_no[k] <- mean(mapply(check_coverage, split(spi_interval, row(spi_interval)), y[test_idx]))
    }


    if (has_unem) {
      ar_base_un <- tryCatch(
        auto.arima(y[train_idx], xreg = unem[train_idx], max.q = 0, D = 0, seasonal = FALSE,
                   allowmean = FALSE, allowdrift = FALSE),
        error = function(e) list(arma = c(1, 0, 0))
      )
      p1_un <- ar_base_un$arma[1]
      fit_un <- tryCatch(
        arima(y[train_idx], order = c(p1_un, 0, 0),
              xreg = X_aug_unem[train_idx, , drop = FALSE],
              include.mean = FALSE),
        error = function(e) NULL
      )
      if (!is.null(fit_un)) {
        yhat_cal_un <- predict(fit_un, n.ahead = length(cal_idx),
                               newxreg = X_aug_unem[cal_idx, , drop = FALSE])$pred
        real_scores_un <- abs(y[cal_idx] - yhat_cal_un)
        synth_scores_un <- as.numeric(abs(hat_y[cal_idx, , drop = FALSE]))
        thr_un <- spi_threshold(real_scores = real_scores_un,
                                synth_scores = synth_scores_un,
                                alpha = alpha)
        yhat_test_un <- predict(fit_un, n.ahead = h,
                                newxreg = X_aug_unem[test_idx, , drop = FALSE])$pred
        spi_interval_un <- cbind(yhat_test_un - thr_un$q_spi, yhat_test_un + thr_un$q_spi)
        len_un[k] <- mean(spi_interval_un[, 2] - spi_interval_un[, 1])
        cov_un[k] <- mean(mapply(check_coverage, split(spi_interval_un, row(spi_interval_un)), y[test_idx]))
      }
    }
  }

  list(
    coverage_without_unem = cov_no,
    length_without_unem = len_no,
    coverage_with_unem = cov_un,
    length_with_unem = len_un
  )
}
