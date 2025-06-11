rm(list = ls())
gc()
source('./utility.R')
df <- read.csv('./meanTopicDistri_mergeCPI100_Unemployment_df19-23_20250301.csv')
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


####################################################################
####################### Prediction Interval ########################
####################################################################

########### PARALLEL ###########
set.seed(1234)
h_c <- seq(8, 15, 1)  # 预测步长
strong_idx <- c(11)  # 强相关变量的索引
p1 <- 2  # AR 阶数
p2 <- 1  # VAR 阶数
rho_c <- seq(0.1, 0.4, length.out = 4)  # 控制变量的相关性
simu_time <- 500  # 模拟次数

all_coverage_result <- matrix(nrow = 0, ncol = length(h_c))
all_len_result <- matrix(nrow = 0, ncol = length(h_c))

library(foreach)
library(doParallel)
library(doRNG)

# 初始化并行集群
cl <- makeCluster(8)  # 使用 8 个核心
registerDoParallel(cl)
clusterExport(cl, c("mvrnorm", "VARX", "%^%"))  # 明确导出 mvrnorm 函数

# 设置随机种子以保证可重复性
clusterSetRNGStream(cl, 1234)

for (d in 1:length(rho_c)) {
  # 使用并行计算处理每个 simu
  results <- foreach(
    simu = 1:simu_time,
    .combine = function(x, y) {
      list(
        coverage = x$coverage + y$coverage,
        len = x$len + y$len,
        success_count = x$success_count + y$success_count
      )
    },
    .packages = c("forecast", "MASS"),  # 添加 MASS 包
    .export = c("generate_arx_data", "check_coverage", 
                "Prediction_BJ", "Bootstrap_powered_ts",
                "X_m", "X_select", "strong_idx", "p1", "p2", "h_c", "rho_c")
  ) %dorng% {
    # 使用 tryCatch 捕获错误，确保程序不会中断
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
        
        # AR模型
        ar_fit <- forecast::Arima(y[obs_idx], order = c(p1, 0, 0), include.mean = FALSE)
        ar_predictions <- forecast::forecast(ar_fit, h = h)
        ar_interval <- cbind(ar_predictions$lower[, 2], ar_predictions$upper[, 2])
        len_simu[k, 1] <- mean(ar_interval[, 2] - ar_interval[, 1])
        coverage_simu[k, 1] <- mean(mapply(check_coverage, split(ar_interval, row(ar_interval)), y[test_idx]))
        
        # Powered+AR+LDA
        powered_lda <- tryCatch({
          Prediction_BJ(X_m, y, tilde_y_m, p1, p2, obs_idx, test_idx, strong_idx, h)
        }, error = function(e) {
          cat("Error in Powered+AR+LDA for h =", h, ":", conditionMessage(e), "\n")
          return(NULL)  # 如果失败返回 NULL
        })
        
        if (!is.null(powered_lda)) {
          powered_interval <- powered_lda$interval
          len_simu[k, 2] <- mean(powered_interval[, 2] - powered_interval[, 1])
          coverage_simu[k, 2] <- mean(mapply(check_coverage, split(powered_interval, row(powered_interval)), y[test_idx]))
        }
        
        # Boot+AR+LDA
        boot_lda <- tryCatch({
          Bootstrap_powered_ts(X_m, y, tilde_y_m, p1, p2, obs_idx, test_idx, strong_idx, h, B = 500)
        }, error = function(e) {
          cat("Error in Boot+AR+LDA for h =", h, ":", conditionMessage(e), "\n")
          return(NULL)  # 如果失败返回 NULL
        })
        
        if (!is.null(boot_lda)) {
          boot_interval <- boot_lda$interval
          len_simu[k, 3] <- mean(boot_interval[, 2] - boot_interval[, 1])
          coverage_simu[k, 3] <- mean(mapply(check_coverage, split(boot_interval, row(boot_interval)), y[test_idx]))
        }
      }
      
      list(coverage = coverage_simu, len = len_simu, success_count = 1)  # 成功计数 +1
    }, error = function(e) {
      cat("Error in simulation", simu, ":", conditionMessage(e), "\n")
      # 如果整个模拟失败，返回零矩阵和 0 成功计数
      list(coverage = matrix(0, nrow = length(h_c), ncol = 3), 
           len = matrix(0, nrow = length(h_c), ncol = 3),
           success_count = 0)
    })
  }
  
  # 提取成功的模拟次数
  successful_simulations <- results$success_count
  
  # 如果成功的模拟次数为 0，则跳过当前 rho_c
  if (successful_simulations == 0) {
    cat("No successful simulations for rho =", rho_c[d], "\n")
    next
  }
  
  # 计算平均值
  coverage_m <- results$coverage / successful_simulations
  len_m <- results$len / successful_simulations
  
  # 合并结果
  all_coverage_result <- rbind(all_coverage_result, t(coverage_m))
  all_len_result <- rbind(all_len_result, t(len_m))
}

# 关闭并行集群
stopCluster(cl)
write.csv(all_coverage_result, './miss_cov_result.csv', row.names = FALSE)
write.csv(all_len_result, './miss_len_result.csv', row.names = FALSE)




