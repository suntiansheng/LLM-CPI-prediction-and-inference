rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
  library(forecast)
  library(parallel)
  library(MASS)
})


source("utility.R")

df <- read.csv("topics_with_pred_score_final.csv", stringsAsFactors = FALSE)
df <- df %>%
  mutate(date = as.Date(date)) %>%
  arrange(date) %>%
  mutate(month = format(date, "%Y-%m"))

df_month <- df %>%
  dplyr::select(month, topic_4, topic_9) %>%
  group_by(month) %>%
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE)), .groups = "drop")

X_m <- as.matrix(df_month[, c("topic_4", "topic_9")])
X_m <- scale(X_m)
X_select <- X_m[, c(1, 2)]

Average_prediction <- function(y, H) {
  pre_c <- rep(0, H)
  for (h in 1:H) {
    pre_c[h] <- mean(y[length(y):(length(y) - h + 1)])
  }
  pre_c
}

check_coverage <- function(pred_intervals, true_values) {
  coverage <- (true_values >= pred_intervals[1]) & (true_values <= pred_intervals[2])
  as.integer(coverage)
}



phi1 <- c(0.5, -0.3)
phi2 <- matrix(c(rep(0.2, 3), rep(-0.2, 3), rep(-0.1, 3)), ncol = 3, byrow = TRUE)
beta1 <- c(0.7, -0.2)
beta2 <- matrix(c(rep(0.1, 3), rep(-0.1, 3)), ncol = 2, byrow = TRUE)

generate_arx_data <- function(X, rho) {
  n <- nrow(X)
  p1 <- length(phi1)
  p2 <- 1
  k <- length(beta1)
  r <- ncol(phi2)

  sigma <- rho * t(matrix(rep(1, r + 1), nrow = 1)) %*% matrix(rep(1, r + 1), nrow = 1)
  diag(sigma) <- 1
  error_c <- mvrnorm(n, mu = rep(0, r + 1), Sigma = sigma)

  y <- numeric(n)
  for (t in (p1 + 1):n) {
    ar_term <- sum(phi1 * y[(t - 1):(t - p1)])
    x_term <- sum(beta1 * X[t, ])
    y[t] <- ar_term + x_term + error_c[t, 1]
  }

  tilde_y <- matrix(0, nrow = n, ncol = r)
  for (t in 2:n) {
    ar_term <- phi2 %*% tilde_y[t - 1, ]
    x_term <- beta2 %*% t(X[t, , drop = FALSE])
    tilde_y[t, ] <- ar_term + x_term + error_c[t, 2:ncol(error_c)]
  }

  list(y = y, tilde_y = tilde_y)
}


h_c <- seq(8, 15, 1)
strong_idx <- c(1, 2)
p1 <- 2
p2 <- 1
rho_c <- seq(0.1, 0.4, length.out = 4)
simu_time <- 500

MSE_result <- matrix(nrow = 0, ncol = length(h_c))
Sign_result <- matrix(nrow = 0, ncol = length(h_c))

for (d in 1:length(rho_c)) {
  results <- parallel::mclapply(seq_len(simu_time), function(simu) {
    tryCatch({
      simu_data <- generate_arx_data(X_select, rho = rho_c[d])
      y <- simu_data$y
      tilde_y_m <- simu_data$tilde_y

      mse_simu <- matrix(0, nrow = length(h_c), ncol = 4)
      sign_simu <- matrix(0, nrow = length(h_c), ncol = 4)

    for (k in 1:length(h_c)) {
      h <- h_c[k]
      test_idx <- (length(y) - h + 1):length(y)
      obs_idx <- setdiff(1:length(y), test_idx)


      ar_fit <- arima(y[obs_idx], order = c(p1, 0, 0), include.mean = FALSE)
      ar_prediction <- predict(ar_fit, n.ahead = h)$pred
      mse_simu[k, 1] <- sqrt(mean((ar_prediction - y[test_idx])^2))
      sign_simu[k, 1] <- mean(sign(ar_prediction) != sign(y[test_idx]))


      mean_prediction <- rep(y[max(obs_idx)], length(test_idx))
      mse_simu[k, 2] <- sqrt(mean((mean_prediction - y[test_idx])^2))
      sign_simu[k, 2] <- mean(sign(mean_prediction) != sign(y[test_idx]))


      ave_prediction <- Average_prediction(y[obs_idx], length(test_idx))
      mse_simu[k, 3] <- sqrt(mean((ave_prediction - y[test_idx])^2))
      sign_simu[k, 3] <- mean(sign(ave_prediction) != sign(y[test_idx]))


      powered_prediction <- LLM_TS.Predict(X_m, y, tilde_y_m, p1, p2,
                                           obs_idx, test_idx, strong_idx, h)
      mse_simu[k, 4] <- sqrt(mean((powered_prediction - y[test_idx])^2))
      sign_simu[k, 4] <- mean(sign(powered_prediction) != sign(y[test_idx]))
    }

      list(mse = mse_simu, sign = sign_simu)
    }, error = function(e) list(error = conditionMessage(e)))
  }, mc.cores = 8, mc.set.seed = TRUE)

  ok_results <- Filter(function(x) is.list(x) && !is.null(x$mse) && !is.null(x$sign), results)
  if (!length(ok_results)) {
    err_msgs <- unique(vapply(results, function(x) if (!is.null(x$error)) x$error else NA_character_, ""))
    err_msgs <- err_msgs[!is.na(err_msgs)]
    stop("All simulations failed for rho = ", rho_c[d], ". Example error: ",
         if (length(err_msgs)) err_msgs[1] else "unknown error",
         call. = FALSE)
  }

  mse_sum <- Reduce(`+`, lapply(ok_results, `[[`, "mse"))
  sign_sum <- Reduce(`+`, lapply(ok_results, `[[`, "sign"))
  denom <- length(ok_results)
  mse_m <- mse_sum / denom
  sign_m <- sign_sum / denom

  MSE_result <- rbind(MSE_result, t(mse_m))
  Sign_result <- rbind(Sign_result, t(sign_m))
}

colnames(Sign_result) <- h_c

method_labels_mse <- rep(c("AR", "RW", "AVE", "LLM-TS"), times = length(rho_c))
rho_labels_mse <- rep(rho_c, each = 4)
rownames(MSE_result) <- paste0("rho_", rho_labels_mse, "_", method_labels_mse)
rownames(Sign_result) <- paste0("rho_", rho_labels_mse, "_", method_labels_mse)

write.csv(MSE_result, "true_MSE.csv", row.names = TRUE)
write.csv(Sign_result, "true_Sign.csv", row.names = TRUE)


all_coverage_result <- matrix(nrow = 0, ncol = length(h_c))
all_len_result <- matrix(nrow = 0, ncol = length(h_c))

for (d in 1:length(rho_c)) {
  results <- parallel::mclapply(seq_len(simu_time), function(simu) {
    tryCatch({
      simu_data <- generate_arx_data(X_select, rho = rho_c[d])
      y <- simu_data$y
      tilde_y_m <- simu_data$tilde_y

      coverage_simu <- matrix(0, nrow = length(h_c), ncol = 3)
      len_simu <- matrix(0, nrow = length(h_c), ncol = 3)

    for (k in 1:length(h_c)) {
      h <- h_c[k]
      test_idx <- (length(y) - h + 1):length(y)
      obs_idx <- setdiff(1:length(y), test_idx)


      ar_fit <- forecast::Arima(y[obs_idx], order = c(p1, 0, 0), include.mean = FALSE)
      ar_predictions_fit <- forecast::forecast(ar_fit, h = length(test_idx))
      ar_interval <- cbind(ar_predictions_fit$lower[, 2], ar_predictions_fit$upper[, 2])
      len_simu[k, 1] <- mean(as.numeric(ar_interval[, 2] - ar_interval[, 1]))
      coverage_simu[k, 1] <- mean(mapply(check_coverage, split(ar_interval, row(ar_interval)), y[test_idx]))


      powered_lda_interval <- LLM_TS.BJ(X_m, y, tilde_y_m, p1, p2,
                                        obs_idx, test_idx, strong_idx, h)$interval
      len_simu[k, 2] <- mean(powered_lda_interval[, 2] - powered_lda_interval[, 1])
      coverage_simu[k, 2] <- mean(mapply(check_coverage, split(powered_lda_interval, row(powered_lda_interval)), y[test_idx]))


      boot_lda_interval <- LLM_TS.BOOT(X_m, y, tilde_y_m, p1, p2,
                                       obs_idx, test_idx, strong_idx, h, B = 500)$interval
      len_simu[k, 3] <- mean(boot_lda_interval[, 2] - boot_lda_interval[, 1])
      coverage_simu[k, 3] <- mean(mapply(check_coverage, split(boot_lda_interval, row(boot_lda_interval)), y[test_idx]))
    }

      list(coverage = coverage_simu, len = len_simu)
    }, error = function(e) list(error = conditionMessage(e)))
  }, mc.cores = 8, mc.set.seed = TRUE)

  ok_results <- Filter(function(x) is.list(x) && !is.null(x$coverage) && !is.null(x$len), results)
  if (!length(ok_results)) {
    err_msgs <- unique(vapply(results, function(x) if (!is.null(x$error)) x$error else NA_character_, ""))
    err_msgs <- err_msgs[!is.na(err_msgs)]
    stop("All simulations failed for rho = ", rho_c[d], ". Example error: ",
         if (length(err_msgs)) err_msgs[1] else "unknown error",
         call. = FALSE)
  }

  coverage_sum <- Reduce(`+`, lapply(ok_results, `[[`, "coverage"))
  len_sum <- Reduce(`+`, lapply(ok_results, `[[`, "len"))
  denom <- length(ok_results)

  all_coverage_result <- rbind(all_coverage_result, t(coverage_sum / denom))
  all_len_result <- rbind(all_len_result, t(len_sum / denom))
}

method_labels_int <- rep(c("AR", "BJ", "BOOT"), times = length(rho_c))
rho_labels_int <- rep(rho_c, each = 3)
rownames(all_coverage_result) <- paste0("rho_", rho_labels_int, "_", method_labels_int)
rownames(all_len_result) <- paste0("rho_", rho_labels_int, "_", method_labels_int)

write.csv(all_coverage_result, "true_Cov.csv", row.names = TRUE)
write.csv(all_len_result, "true_Len.csv", row.names = TRUE)
