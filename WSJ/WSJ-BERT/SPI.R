













library(MTS)
library(stats)

quiet_eval <- function(expr) {
    invisible(capture.output(result <- eval(expr, envir = parent.frame())))
    result
}




spi_transport_scores <- function(real_scores, synth_scores) {
    u <- ecdf(real_scores)(real_scores)
    quantile(synth_scores, probs = u, type = 8, names = FALSE)
}




spi_threshold <- function(real_scores, synth_scores, alpha = 0.05) {
    z_real  <- spi_transport_scores(real_scores, synth_scores)


    q_synth <- quantile(c(synth_scores, z_real), probs = 1 - alpha, type = 8)


    idx <- which(z_real <= q_synth)
    if (length(idx) == 0) {
        q_spi <- max(real_scores)
    } else {
        q_spi <- max(real_scores[idx])
    }

    list(q_spi = q_spi,
         q_synth = q_synth,
         z_real  = z_real)
}




recursive_arx_forecast <- function(y, X, coef_ar, coef_x, p, start_idx, h) {

    if (start_idx + h > nrow(X)) {
        stop("Not enough future surrogate rows in X for recursive forecast.")
    }
    y_ext <- y[1:start_idx]
    yhat  <- numeric(h)

    for (k in 1:h) {


        x_future <- X[start_idx + k, ]


        lags <- rev(tail(y_ext, p))

        ar_term <- sum(coef_ar * lags)
        x_term  <- sum(coef_x * x_future)

        yhat[k] <- ar_term + x_term


        y_ext <- c(y_ext, yhat[k])
    }
    yhat
}




spi_recursive_one_h <- function(
        y, X, tilde_y, X_sur, strong_idx,
        obs_idx, test_idx, h, p2 = 1, alpha = 0.05) {




    y_train <- y[obs_idx]
    X_train <- X[obs_idx, strong_idx, drop = FALSE]

    ar_base <- tryCatch(
        quiet_eval(quote(auto.arima(y_train, max.q = 0, seasonal = FALSE,
                   allowmean = FALSE, allowdrift = FALSE))),
        error = function(e) list(arma = c(1,0,0))
    )
    p_ar <- ar_base$arma[1]

    fit_ar <- quiet_eval(quote(arima(
        y_train,
        order = c(p_ar, 0, 0),
        xreg  = X_train,
        include.mean = FALSE
    )))

    coef <- fit_ar$coef
    coef_ar <- coef[1:p_ar]
    coef_x  <- coef[(p_ar+1):length(coef)]




    yhat_1step <- fitted(fit_ar)
    real_scores <- abs(y_train - yhat_1step)




    tilde_train <- tilde_y[obs_idx, , drop = FALSE]
    X_train_sur <- X_sur[obs_idx, strong_idx, drop = FALSE]

    varx_fit <- quiet_eval(quote(VARX(
        zt = tilde_train,
        p  = p2,
        xt = X_train_sur,
        include.mean = FALSE
    )))


    Phi_arr <- varx_fit$Phi
    Phi_list <- if (length(dim(Phi_arr)) == 3) {
        lapply(seq_len(dim(Phi_arr)[3]), function(j) as.matrix(Phi_arr[,,j]))
    } else {
        list(as.matrix(Phi_arr))
    }


    synth_scores <- rowMeans(abs(varx_fit$residuals))




    q_spi <- spi_threshold(real_scores, synth_scores, alpha)$q_spi




    start_idx <- tail(obs_idx, 1)
    yhat_h <- recursive_arx_forecast(
        y = y,
        X = X[, strong_idx, drop = FALSE],
        coef_ar = coef_ar,
        coef_x  = coef_x,
        p = p_ar,
        start_idx = start_idx,
        h = h
    )

    list(
        pred = yhat_h,
        lower = yhat_h - q_spi,
        upper = yhat_h + q_spi,
        q_spi = q_spi
    )
}




spi_recursive_multi <- function(
        y, X, tilde_y, strong_idx,
        horizon_list = 1:12, p2 = 1, alpha = 0.05) {

    n <- length(y)
    max_h <- max(horizon_list)
    if (max_h >= n) {
        stop("Max horizon must be smaller than series length.")
    }
    X_sur <- X

    out <- list()
    for (h in horizon_list) {
        obs_idx <- seq_len(n - h)
        test_idx <- (n - h + 1):n

        spi_out <- spi_recursive_one_h(
            y = y,
            X = X,
            tilde_y = tilde_y,
            X_sur = X_sur,
            strong_idx = strong_idx,
            obs_idx = obs_idx,
            test_idx = test_idx,
            h = h,
            p2 = p2,
            alpha = alpha
        )
        truth <- y[test_idx]
        coverage <- mean(truth >= spi_out$lower & truth <= spi_out$upper)

        out[[paste0("h", h)]] <- list(
            pred = spi_out$pred,
            lower = spi_out$lower,
            upper = spi_out$upper,
            q_spi = spi_out$q_spi,
            rmse  = sqrt(mean((spi_out$pred - truth)^2)),
            sign_acc = mean(sign(spi_out$pred) == sign(truth)),
            coverage = coverage,
            length   = mean(spi_out$upper - spi_out$lower)
        )
    }

    out
}




spi_multi_horizon <- function(y, X_m, tilde_y_m, strong_idx, h_vec,
                              unem, has_unem = FALSE,
                              cal_frac = 0.2, alpha = 0.05, p2 = 1) {
    if (length(h_vec) == 0) {
        return(list(
            coverage_without_unem = numeric(),
            length_without_unem = numeric(),
            coverage_with_unem = numeric(),
            length_with_unem = numeric()
        ))
    }

    n <- length(y)
    res_no <- spi_recursive_multi(
        y = y,
        X = X_m,
        tilde_y = tilde_y_m,
        strong_idx = strong_idx,
        horizon_list = h_vec,
        p2 = p2,
        alpha = alpha
    )

    coverage_without_unem <- numeric(length(h_vec))
    length_without_unem <- numeric(length(h_vec))
    for (i in seq_along(h_vec)) {
        out <- res_no[[paste0("h", h_vec[i])]]
        truth <- y[(n - h_vec[i] + 1):n]
        coverage_without_unem[i] <- mean(truth >= out$lower & truth <= out$upper, na.rm = TRUE)
        length_without_unem[i] <- mean(out$upper - out$lower, na.rm = TRUE)
    }

    coverage_with_unem <- rep(NA_real_, length(h_vec))
    length_with_unem <- rep(NA_real_, length(h_vec))
    if (has_unem) {
        X_with_unem <- cbind(X_m, unem)
        strong_idx_unem <- unique(c(strong_idx, ncol(X_with_unem)))
        res_un <- spi_recursive_multi(
            y = y,
            X = X_with_unem,
            tilde_y = tilde_y_m,
            strong_idx = strong_idx_unem,
            horizon_list = h_vec,
            p2 = p2,
            alpha = alpha
        )
        for (i in seq_along(h_vec)) {
            out <- res_un[[paste0("h", h_vec[i])]]
            truth <- y[(n - h_vec[i] + 1):n]
            coverage_with_unem[i] <- mean(truth >= out$lower & truth <= out$upper, na.rm = TRUE)
            length_with_unem[i] <- mean(out$upper - out$lower, na.rm = TRUE)
        }
    }

    list(
        coverage_without_unem = coverage_without_unem,
        length_without_unem = length_without_unem,
        coverage_with_unem = coverage_with_unem,
        length_with_unem = length_with_unem
    )
}
