# Load library
library(dplyr)
library(evd)         # fitting GPD
library(POT)         # Peaks Over Threshold
library(ismev)       # MLE GPD
library(ggplot2)     # visualisation
library(tidyr)
library(glmnet)      # fused ridge regression
library(VineCopula)  # copula modeling
library(copula)
library(numDeriv)    # matriks Hessian
library(MASS)        # ginv()
library(purrr)
library(igraph)
library(evir)
library(ggspatial)
library(sf)
library(akima)       # interpolasi IDW
library(splines)
library(viridis)
library(goftest)
library(geosphere)   # distm/jarak euclidian
library(reshape2)    # Gaussian Copula
library(tseries)     # correlation rainfall 
library(fExtremes)   # extremal index
library(maps)        # maps
library(gstat)
library(stars)
library(stringr)




# =========================
#          DATA
# =========================
data <- read.csv("Data RAINFALL Sumut 14-24 Complete.csv")  
data <- data %>%
  filter(Station %in% c("FL Tobing", "Aek Gondang", "Kualanamu", "Maritim Belawan", "Binaka"))
head(data)



# ====================
#   DESCRIPTION DATA
# ====================
# Statistik deskriptif
summary_data <- data %>%
  group_by(Station) %>%  
  summarise(
    Min = min(Rainfall, na.rm = TRUE),
    Q1 = quantile(Rainfall, 0.25, na.rm = TRUE),
    Median = median(Rainfall, na.rm = TRUE),
    Mean = mean(Rainfall, na.rm = TRUE),
    Q3 = quantile(Rainfall, 0.75, na.rm = TRUE),
    Max = max(Rainfall, na.rm = TRUE),
    SD  = sd(Rainfall, na.rm = TRUE)
  )
print(summary_data)




# ========================
#   THRESHOLD SELECTION
# ========================
estimasi_gpd_location <- function(df, location) {
  data_rain <- df %>% filter(Station == location) %>% pull(Rainfall)
  percentiles <- seq(0.50, 0.95, by = 0.01)
  thresholds <- quantile(data_rain, percentiles, na.rm = TRUE)
  
  results <- data.frame()
  
  for (i in seq_along(thresholds)) {
    u <- thresholds[i]
    if (sum(data_rain > u) < 10) next
    
    fit_gpd <- ismev::gpd.fit(data_rain, threshold = u, show = FALSE)
    if (is.null(fit_gpd)) next
    
    data_extreme <- fit_gpd$data[fit_gpd$data > u] - u
    
    xi <- fit_gpd$mle[2]
    sigma <- fit_gpd$mle[1]
    data_norm <- data_extreme / sigma
    
    if (is.na(sigma) || sigma <= 0 || is.na(xi)) next
    if (any(data_norm < 0)) next
    
    ad <- tryCatch({
      goftest::ad.test(data_norm, null = function(x) {
        inside <- 1 + xi * x
        if (any(inside <= 0)) return(rep(NA, length(x)))
        return(1 - inside^(-1 / xi))
      })
    }, error = function(e) return(NULL))
    
    if (is.null(ad) || is.na(ad$statistic) || is.na(ad$p.value)) next
    
    n_exc <- length(data_extreme)
    
    results <- rbind(results, data.frame(
      Percentile = percentiles[i] * 100,
      Threshold = round(u, 2),
      xi = round(xi, 2),
      sigma = round(sigma, 2),
      AD_Stat = round(ad$statistic, 3),
      AD_p_value = round(ad$p.value, 3),
      Total_Exceedance = n_exc,
      Station = location
    ))
  }
  
  results <- results %>% arrange(Threshold)
  results <- results %>% mutate(
    MRL = map_dbl(Threshold, function(u) {
      exc <- data_rain[data_rain > u]
      if (length(exc) > 0) mean(exc - u) else NA_real_
    }),
    MRL_diff = c(NA, diff(MRL)),
    MRL_stabilitas = abs(MRL_diff)
  )
  
  results_filtered <- results %>%
    filter(between(xi, -0.5, 0.5), !is.na(AD_p_value))
  
  best <- results_filtered %>%
    arrange(desc(AD_p_value), MRL_stabilitas) %>%
    slice(1) %>%
    mutate(Threshold_Best = TRUE)
  
  results <- results %>%
    mutate(Threshold_Best = Threshold %in% best$Threshold)
  
  return(results)
}

# MRL Plot untuk semua Stasiun
mrl_all <- map_df(unique(data$Station), function(station_name) {
  data_station <- data %>% filter(Station == station_name)
  
  th_seq <- quantile(data_station$Rainfall,
                     probs = seq(0.50, 0.95, 0.01),
                     na.rm = TRUE)
  
  mrl <- map_dbl(th_seq, function(u) {
    exc <- data_station$Rainfall[data_station$Rainfall > u]
    if (length(exc) > 0) mean(exc - u) else NA_real_
  })
  
  data.frame(
    Station = station_name,
    Threshold = th_seq,
    MRL = mrl
  )
})

threshold_summary <- unique(data$Station) %>%
  map_df(~estimasi_gpd_location (data, .x)) %>%
  dplyr::filter(Threshold_Best)

mrl_all <- mrl_all %>%
  left_join(
    threshold_summary %>% dplyr::select(Station, Threshold),
    by = "Station",
    suffix = c("", "_final")
  )

ggplot(mrl_all, aes(x = Threshold, y = MRL)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(color = "steelblue", size = 1.5) +
  
  geom_vline(aes(xintercept = Threshold_final),
             linetype = "dashed",
             color = "red",
             linewidth = 0.8) +
  facet_wrap(~Station, scales = "free") +
  labs(
    title = "Mean Residual Life (MRL) Plot",
    x = "Threshold (mm)",
    y = "Mean Excess"
  ) +
  theme_minimal(base_size = 13)

print(threshold_summary)




# ====================
#   DATA EXCEEDANCE
# ====================
data_exceedance <- data %>%
  left_join(threshold_summary %>% 
              rename(Threshold = Threshold) %>% 
              dplyr::select(Station, Threshold), 
            by = "Station") %>%
  filter(Rainfall > Threshold) %>%  
  mutate(
    Excess = Rainfall - Threshold
  ) %>%
  dplyr::select(Station, Rainfall, Threshold, Excess, Latitude, Longitude)
head(data_exceedance)

data_exceedance %>%
  group_by(Station) %>%
  summarise(n_exceed = n(), .groups = "drop")




# ===============================
#   LOG-LIKELIHOOD GPD FUNCTION 
# ===============================
# Shape  Parameter = konstan
xi_data <- data_exceedance$Excess
fit_xi <- ismev::gpd.fit(xi_data, threshold = 0, show = FALSE)
xi_konstan <- fit_xi$mle[2]
print(xi_konstan)

data_list <- split(data_exceedance$Excess,
                   data_exceedance$Station)
loglik_gpd <- function(sigma, xi, data_list){
  total_loglik <- 0
  
  for(i in seq_along(data_list)){
    y <- data_list[[i]]
    n <- length(y)
    
    loglik_i <-
      -n*log(sigma[i]) -
      (1/xi + 1) *
      sum(log(1 + xi*y/sigma[i]))
    
    total_loglik <- total_loglik + loglik_i
  }
  
  return(total_loglik)
}




# ==================
#   FUSED PENALTY
# ==================
location_coords <- data_exceedance %>%
  group_by(Station) %>%
  summarise(
    Latitude = mean(Latitude, na.rm = TRUE),
    Longitude = mean(Longitude, na.rm = TRUE)
  )

coords_matrix <- as.matrix(location_coords[, c("Longitude", "Latitude")])
long_matrix <- distm(coords_matrix, fun = distHaversine)

g <- graph_from_adjacency_matrix(
  long_matrix,
  mode = "undirected",
  weighted = TRUE,
  diag = FALSE
)
V(g)$name <- location_coords$Station            
mst_g <- mst(g, weights = E(g)$weight)          
edge_df <- as_data_frame(mst_g, what = "edges") 

Wpq_matrix <- 1 / long_matrix
diag(Wpq_matrix) <- 0  
Wpq_matrix_scaled <- Wpq_matrix / mean(Wpq_matrix[Wpq_matrix > 0])
edge_df$Wpq <- mapply(
  function(a, b){
    i <- match(a, location_coords$Station)
    j <- match(b, location_coords$Station)
    
    Wpq_matrix_scaled[i, j]
  },
  edge_df$from,
  edge_df$to
)
edge_df <- edge_df %>%
  rename(
    Location_1 = from,
    Location_2 = to
  )
print(edge_df)

data_list <- data_exceedance %>%
  group_by(Station) %>%
  summarise(Excess = list(Excess), .groups = "drop")

input_fused_penalty <- list(
  xi = xi_konstan,
  location = location_coords,
  edge = edge_df,
  W_matrix = Wpq_matrix_scaled,
  excess_location = data_list
)

saveRDS(input_fused_penalty, "input_fused_penalty.rds")  

# LOAD DATA INPUT
input <- readRDS("input_fused_penalty.rds")

data_list <- input$excess_location$Excess
n_lokasi <- length(data_list)
station_names <- input$location$Station
xi_konstan <- input$xi
edge_df <- input$edge

estimation_gpd_without_regularization <- function(df, location){
  data_excess <- df %>%
    filter(Station == location)
  
  excess <- data_excess$Excess
  
  fit <- ismev::gpd.fit(
    excess,
    threshold = 0,
    show = FALSE
  )
  
  data.frame(
    Station = location,
    Sigma = fit$mle[1],
    Xi = fit$mle[2]
  )
}
station_list <- unique(data_exceedance$Station)
result_estimation <- map_dfr(station_list, ~ estimation_gpd_without_regularization(data_exceedance, .x))
print(result_estimation)
sigma_init <- input$excess_location %>%
  left_join(result_estimation, by = "Station") %>%
  pull(Sigma)
identical(input$excess_location$Station, 
          result_estimation$Stasiun[match(input$excess_location$Station, result_estimation$Station)])   

print(input_fused_penalty)


# ===============
#   FUSED LASSO
# ===============
penalty_lasso <- function(sigma, edge_df, station_names, lambda){
  penalty <- 0
  
  for(k in 1:nrow(edge_df)){
    p <- match(edge_df$Location_1[k], station_names)
    q <- match(edge_df$Location_2[k], station_names)
    
    penalty <- penalty +
      edge_df$Wpq[k] *
      abs(sigma[p] - sigma[q])
  }
  
  return(lambda * penalty)
}

objective_lasso <- function(sigma,
                            data_list,
                            xi,
                            edge_df,
                            station_names,
                            lambda){
  
  ll <- loglik_gpd(sigma, xi, data_list)
  
  pen_lasso <- penalty_lasso(
    sigma,
    edge_df,
    station_names,
    lambda
  )
  
  return(-(ll - pen_lasso))
}

lambda_grid <- seq(0.01, 20, by = 0.5)
tic_lambda_lasso <- data.frame(
  lambda = lambda_grid,
  TIC_lambda = NA
)

for(i in seq_along(lambda_grid)){
  lambda <- lambda_grid[i]
  
  res <- optim(
    par = sigma_init,
    fn = objective_lasso,
    method = "L-BFGS-B",
    lower = rep(1e-6, n_location),
    data_list = data_list,
    xi = xi_konstan,
    edge_df = edge_df,
    station_names = station_names,
    lambda = lambda
  )
  
  sigma_tmp <- res$par
  loglik_tmp <- -objective_lasso(
    sigma_tmp,
    data_list,
    xi_konstan,
    edge_df,
    station_names,
    lambda
  )
  hess_tmp <- hessian(
    func = objective_lasso,
    x = sigma_tmp,
    data_list = data_list,
    xi = xi_konstan,
    edge_df = edge_df,
    station_names = station_names,
    lambda = lambda
  )
  grad_tmp <- grad(
    func = objective_lasso,
    x = sigma_tmp,
    data_list = data_list,
    xi = xi_konstan,
    edge_df = edge_df,
    station_names = station_names,
    lambda = lambda
  )
  H_inv <- tryCatch(solve(hess_tmp), error = function(e) MASS::ginv(hess_tmp))
  J_matrix <- grad_tmp %*% t(grad_tmp)
  
  tic_lambda_lasso$TIC_lambda[i] <- -2 * (loglik_tmp - sum(diag(H_inv %*% J_matrix)))
}

best_lambda_lasso <- tic_lambda_lasso$lambda[which.min(tic_lambda_lasso$TIC_lambda)]
print(best_lambda_lasso)

res_best_lasso <- optim(
  par = sigma_init,
  fn = objective_lasso,
  method = "L-BFGS-B",
  lower = rep(1e-6, n_location),
  data_list = data_list,
  xi = xi_konstan,
  edge_df = edge_df,
  station_names = station_names,
  lambda = best_lambda_lasso,
  control = list(maxit = 5000, factr = 1e-10)
)

sigma_hat_lasso <- res_best_lasso$par
result_lasso <- data.frame(
  Station = station_names,
  Sigma_FusedLasso = sigma_hat_lasso
)
print(result_lasso)


# ===============
#   FUSED RIDGE
# ===============
penalty_ridge <- function(sigma, edge_df, station_names, lambda){
  penalty <- 0
  
  for(k in 1:nrow(edge_df)){
    p <- match(edge_df$Location_1[k], station_names)
    q <- match(edge_df$Location_2[k], station_names)
    
    penalty <- penalty +
      edge_df$Wpq[k] *
      (sigma[p] - sigma[q])^2
  }
  
  return(lambda * penalty)
}

# Fungsi objective ridge untuk optim
objective_ridge <- function(sigma,
                            data_list,
                            xi,
                            edge_df,
                            station_names,
                            lambda){
  
  ll <- loglik_gpd(sigma, xi, data_list)
  
  pen_ridge <- penalty_ridge(
    sigma,
    edge_df,
    station_names,
    lambda
  )
  
  return(-(ll - pen_ridge))
}

lambda_grid <- seq(0.01, 20, by = 0.5)
tic_lambda_ridge <- data.frame(
  lambda = lambda_grid,
  TIC_lambda = NA
)

for(i in seq_along(lambda_grid)){
  lambda <- lambda_grid[i]
  
  res <- optim(
    par = sigma_init,
    fn = objective_ridge,
    method = "L-BFGS-B",
    lower = rep(1e-6, n_lokasi),
    data_list = data_list,
    xi = xi_konstan,
    edge_df = edge_df,
    station_names = station_names,
    lambda = lambda
  )
  
  sigma_tmp <- res$par
  loglik_tmp <- -objective_ridge(
    sigma_tmp,
    data_list,
    xi_konstan,
    edge_df,
    station_names,
    lambda
  )
  hess_tmp <- hessian(
    func = objective_ridge,
    x = sigma_tmp,
    data_list = data_list,
    xi = xi_konstan,
    edge_df = edge_df,
    station_names = station_names,
    lambda = lambda
  )
  grad_tmp <- grad(
    func = objective_ridge,
    x = sigma_tmp,
    data_list = data_list,
    xi = xi_konstan,
    edge_df = edge_df,
    station_names = station_names,
    lambda = lambda
  )
  H_inv <- tryCatch(solve(hess_tmp), error = function(e) MASS::ginv(hess_tmp))
  J_matrix <- grad_tmp %*% t(grad_tmp)
  
  tic_lambda_ridge$TIC_lambda[i] <- -2 * (loglik_tmp - sum(diag(H_inv %*% J_matrix)))
}

best_lambda_ridge <- tic_lambda_ridge$lambda[which.min(tic_lambda_ridge$TIC_lambda)]
print(best_lambda_ridge)

res_best_ridge <- optim(
  par = sigma_init,
  fn = objective_ridge,
  method = "L-BFGS-B",
  lower = rep(1e-6, n_location),
  data_list = data_list,
  xi = xi_konstan,
  edge_df = edge_df,
  station_names = station_names,
  lambda = best_lambda_ridge,
  control = list(maxit = 5000, factr = 1e-10)
)

sigma_hat_ridge <- res_best_ridge$par
result_ridge <- data.frame(
  Station = station_names,
  Sigma_FusedRidge = sigma_hat_ridge
)
print(result_ridge)



                    
# ============================================
#   PARAMETER GAUSSIAN COPULA (R) ESTIMATION
# ============================================
fix_station <- function(x) {
  x %>%
    str_replace_all(" ", "\\.") %>%  
    stringr::str_squish()           
}
threshold_summary <- threshold_summary %>%
  dplyr::mutate(Station = fix_station(Station))
result_test_gpd <- result_test_gpd %>%
  dplyr::mutate(Station = fix_station(Station))
data_exceedance <- data_exceedance %>%
  dplyr::mutate(Station = fix_station(Station))

transform_uniform_semiparametric <- function(data_wide,
                                             parameter_gpd,
                                             threshold_summary){
  result <- data_wide
  name_station <- colnames(data_wide)[-1]
  
  for(st in name_station){
    x <- data_wide[[st]]
    threshold <- threshold_summary %>%
      filter(Station == st) %>%
      pull(Threshold)
    sigma <- parameter_gpd %>%
      filter(Station == st) %>%
      pull(Sigma_Hat)
    xi <- parameter_gpd %>%
      filter(Station == st) %>%
      pull(Xi_Konstan)
    
    lambda_u <- mean(x > threshold)
    
    # empirical CDF
    F_emp <- ecdf(x)
    u <- numeric(length(x))
    
    for(i in seq_along(x)){
      if(x[i] <= threshold){
        u[i] <- (1 - lambda_u) * F_emp(x[i])
      } else {
        y <- x[i] - threshold
        gpd_tail <- evd::pgpd(
          q = y,
          loc = 0,
          scale = sigma,
          shape = xi
        )
        u[i] <- (1 - lambda_u) + lambda_u * gpd_tail
      }
    }
    
    result[[st]] <- u
  }
  
  return(hasil)
}

parameter_lasso <- result_test_gpd %>%
  filter(Method == "Fused Lasso")
data_uniform_lasso <- transform_uniform_semiparametric(
  data_wide,
  parameter_lasso,
  threshold_summary)

parameter_ridge <- result_test_gpd %>%
  filter(Method == "Fused Ridge")
data_uniform_ridge <- transform_uniform_semiparametric(
  data_wide,
  parameter_ridge,
  threshold_summary)

# Gaussian Copula Function
gauss_copula <- function(data_uniform, name_method){
  u_matrix <- as.matrix(data_uniform[,-1])
  
  z_matrix <- qnorm(u_matrix)
  
  R_hat <- cor(z_matrix)
  
  cop_gauss <- normalCopula(
    param = P2p(R_hat),
    dim = ncol(R_hat),
    dispstr = "un"
  )
  
  return(list(
    R_hat = R_hat,
    copula_model = cop_gauss
  ))
}

copula_lasso <- gauss_copula(data_uniform_lasso, "Fused Lasso")
copula_ridge <- gauss_copula(data_uniform_ridge, "Fused Ridge")

R_hat_lasso <- copula_lasso$R_hat
R_hat_ridge <- copula_ridge$R_hat

cop_model_lasso <- copula_lasso$copula_model
cop_model_ridge <- copula_ridge$copula_model

print(list(R_hat_lasso, R_hat_ridge))
print(list(cop_model_lasso, cop_model_ridge))

saveRDS(copula_lasso, "copula_lasso.rds")
saveRDS(copula_ridge, "copula_ridge.rds")




# ================
#   RETURN LEVEL
# ================
count_return_level <- function(exceedance, 
                               parameter_gpd,
                               threshold,
                               T_year){
  
  T_day <- T_year * 365
  
  lambda_table <- exceedance %>%
    group_by(Station) %>%
    summarise(lambda_u = n()/4018, .groups = "drop")
  
  param <- parameter_gpd %>%
    dplyr::select(Station, Sigma_Hat, Xi_Konstan) %>%
    left_join(threshold %>%
                dplyr::select(Station, Threshold),
              by = "Station") %>%
    left_join(lambda_table, by="Station") 
  
  probabilitas <- 1 - 1/(T_day * param$lambda_u)
  
  return_uni <- param %>%
    mutate(
      Return_Level = mapply(
        function(p, sigma, xi, threshold){
          if(is.na(xi) || is.na(sigma) || is.na(threshold)){
            return(NA)}
          if(abs(xi) < 1e-6){
            threshold - sigma * log(1 - p)
          } else {
            threshold + (sigma / xi) * ((1 - p)^(-xi) - 1)}
        },
        p = probabilitas,
        sigma = Sigma_Hat,
        xi = Xi_Konstan,
        threshold = Threshold
      )
    )
  
  return(return_uni)
}

period_list <- c(10, 25, 100)
result_RL <- list()

for (T in period_list) {
  ReturnLvl_Lasso <- count_return_level(data_exceedance, 
                                        parameter_lasso,
                                        threshold_summary,
                                        T_year = T)
  ReturnLvl_Ridge <- count_return_level(data_exceedance,
                                        parameter_ridge,
                                        threshold_summary,
                                        T_year = T)
  
  df_temp <- ReturnLvl_Lasso %>%
    dplyr::select(Station, RL_Lasso = Return_Level) %>%
    left_join(
      ReturnLvlRidge %>%
        dplyr::select(Station, RL_Ridge = Return_Level),
      by = "Station"
    )
  df_temp$Period_Year <- T
  
  result_RL[[as.character(T)]] <- df_temp
}

return_level_final <- bind_rows(result_RL)
print(return_level_final)




# ==========================================
#   COMPARASION METHOD FUSED PENLATY (TIC)
# ==========================================
loglik_gpd <- function(x, sigma, xi) {
  if (any(sigma <= 0)) return(-Inf)
  
  cond <- 1 + xi * x / sigma
  if (any(cond <= 0)) return(-Inf)
  if (abs(xi) < 1e-6) {
    ll <- -sum(log(sigma) + x / sigma)
  } else {
    ll <- sum(-log(sigma) - (1/xi + 1) * log(cond))
  }
  return(ll)
}

loglik_gpd_i <- function(x, sigma, xi) {
  if (sigma <= 0) return(-Inf)
  cond <- 1 + xi * x / sigma
  if (cond <= 0) return(-Inf)
  if (abs(xi) < 1e-6) {
    ll <- -(log(sigma) + x / sigma)
  } else {
    ll <- -(log(sigma) + (1/xi + 1) * log(cond))
  }
  return(ll)
}

count_TIC_Uni <- function(data, sigma_vector, xi) {
  station_unique <- unique(data$Station)
  n_station <- length(sigma_vector)
  
  data$Station_index <- match(data$Station, station_unique)
  
  sigma_obs <- sigma_vector[data$Station_index]
  
  loglik <- loglik_gpd(
    x = data$Excess,
    sigma = sigma_obs,
    xi = xi
  )
  
  param <- c(sigma_vector, xi)
  
  grad_func <- function(p){
    sigma_temp <- p[1:n_stasiun]
    xi_temp <- p[length(p)]
    sigma_obs_temp <- sigma_temp[data$Station_index]
    
    return(loglik_gpd(x = data$Excess, sigma = sigma_obs_temp, xi = xi_temp))
  }
  
  J <- -hessian(grad_func, param)
  
  # Hitung Score Matrix (K)
  score_matrix <- matrix(0, nrow = nrow(data), ncol = length(param))
  
  for(i in 1:nrow(data)){
    idx_station <- data$Station_index[i]
    
    score_i <- grad(
      function(p){
        sigma_temp <- p[1:n_station]
        xi_temp <- p[length(p)]
        
        sigma_i <- sigma_temp[idx_stasiun]
        
        return(loglik_gpd_i(x = data$Excess[i], sigma = sigma_i, xi = xi_temp))
      },
      param
    )
    
    score_matrix[i, ] <- score_i
  }
  
  K <- t(score_matrix) %*% score_matrix
  
  inv_J <- tryCatch(solve(J), error = function(e) {
    warning("Matriks Hessian (J) cannot be directly inverted. Use the pseudo-inverse.")
    return(MASS::ginv(J))
  })
  
  TIC <- -2 * loglik + 2 * sum(diag(inv_J %*% K))
  
  return(list(
    loglik = loglik,
    TIC = TIC
  ))
}
data_exceedance <- data_exceedance %>%
  mutate(Station = str_replace_all(Station, "\\.", " "))
str(data_exceedance)
str(result_estimation)
str(result_lasso)
str(result_ridge)

# Original
data_tic_ori <- data_exceedance %>%
  left_join(result_estimation %>% 
              dplyr::select(Station, Sigma), by = "Station")
sigma_unique_ori <- data_tic_ori %>%
  distinct(Station, .keep_all = TRUE) %>%
  pull(Sigma)
tic_ori <- count_TIC_Uni(
  data = data_tic_ori,
  sigma_vector = sigma_unique_ori,
  xi = unique(result_test_gpd$Xi_Konstan)
)

# Fused Lasso 
data_tic_lasso <- data_exceedance %>%
  left_join(result_lasso %>% 
              dplyr::select(Station, Sigma_FusedLasso), by = "Station")
sigma_unique_lasso <- data_tic_lasso %>%
  distinct(Station, .keep_all = TRUE) %>%
  pull(Sigma_FusedLasso)
tic_lasso <- count_TIC_Uni(
  data = data_tic_lasso,
  sigma_vector = sigma_unique_lasso,
  xi = unique(result_test_gpd$Xi_Konstan)
)

# Fused Ridge 
data_tic_ridge <- data_exceedance %>%
  left_join(result_ridge %>% 
              dplyr::select(Station, Sigma_FusedRidge), by = "Station")
sigma_unique_ridge <- data_tic_ridge %>%
  distinct(Station, .keep_all = TRUE) %>%
  pull(Sigma_FusedRidge)
tic_ridge <- count_TIC_Uni(
  data = data_tic_ridge,
  sigma_vector = sigma_unique_ridge,
  xi = unique(result_test_gpd$Xi_Konstan)
)

comparasion_TIC_Uni <- data.frame(
  Method = c("Fused Lasso", "Fused Ridge", "Original"),
  LogLik = c(tic_lasso$loglik, tic_ridge$loglik, tic_ori$loglik),
  TIC = c(tic_lasso$TIC, tic_ridge$TIC, tic_ori$TIC)
)
print(comparasion_TIC_Uni)




# =======================
#   VISUALISASI HEATMAP
# =======================
data_coordinate <- data %>%
  dplyr::select(Station, Latitude, Longitude) %>%
  distinct()
data_spasial <- return_level_final %>%
  dplyr::mutate(Station = str_replace_all(Station, "\\.", " ")) %>%
  left_join(data_coordinate, by = "Station")

data_long <- data_spasial %>%
  dplyr::select(Station, Longitude, Latitude, Period_Year, RL_Lasso, RL_Ridge) %>%
  pivot_longer(
    cols = c(RL_Lasso, RL_Ridge),
    names_to = "Method",
    values_to = "return_level"
  ) %>%
  mutate(
    Method = recode(
      Method,
      RL_Lasso = "Fused Lasso",
      RL_Ridge = "Fused Ridge"
    ),
    Period_Year = factor(
      Period_Year,
      levels = c(10, 25, 100),
      labels = c("10 Year", "25 Year", "100 Year")
    )
  )

indo <- st_read("DataMaps.shp")
indo <- st_make_valid(indo)
region_select <- indo %>%
  filter(KAB_KOTA %in% c("TAPANULI TENGAH",
                         "PADANG LAWAS UTARA",
                         "DELI SERDANG",
                         "KOTA MEDAN",
                         "NIAS"))
grid <- st_make_grid(region_select, cellsize = 0.05) %>%
  st_as_sf()

interpolasi_idw <- function(df){
  
  df_clean <- df %>%
    dplyr::select(Longitude, Latitude, return_level)
  station_sf <- st_as_sf(df_clean, coords = c("Longitude", "Latitude"), crs = 4326)
  
  idw_result <- gstat::idw(
    return_level ~ 1,
    locations = station_sf,
    newdata = grid
  )
  idw_sf <- st_as_sf(idw_result)

  idw_clip <- st_intersection(idw_sf, st_make_valid(region_select))
  return(idw_clip)
}

idw_all <- data_long %>%
  filter(!is.na(Latitude) & !is.na(Longitude)) %>%
  group_by(Method, Period_Year) %>%
  group_split() %>%
  map_df(~{
    result <- interpolasi_idw(.x)
    result$Method <- unique(.x$Method)
    result$Period_Year <- unique(.x$Period_Year)
    return(result)
  })
idw_all <- st_make_valid(idw_all)

station_sf <- st_as_sf(data_coordinate,
                       coords = c("Longitude", "Latitude"),
                       crs = 4326)

bbox_sumut <- st_bbox(stasiun_sf)
bbox_sumut["xmin"] <- bbox_sumut["xmin"] - 1
bbox_sumut["xmax"] <- bbox_sumut["xmax"] + 1
bbox_sumut["ymin"] <- bbox_sumut["ymin"] - 1
bbox_sumut["ymax"] <- bbox_sumut["ymax"] + 1
indo_sumut <- st_crop(indo, bbox_sumut)
indo_sumut <- st_simplify(indo_sumut, dTolerance = 0.01)

idw_ridge <- idw_all %>%
  filter(Method == "Fused Ridge") %>%
  st_cast("POINT") %>%
  mutate(
    x = st_coordinates(.)[,1],
    y = st_coordinates(.)[,2]
  )
idw_ridge %>%
  group_by(Period_Year) %>%
  summarise(min = min(var1.pred),
            max = max(var1.pred))
idwr_10 <- idw_ridge %>% filter(Period_Year == "10 Year")
idwr_25 <- idw_ridge %>% filter(Period_Year == "25 Year")
idwr_100 <- idw_ridge %>% filter(Period_Year == "100 Year")
st_bbox(indo_sumut)

plot_rl <- function(data_plot, title){
  ggplot() +
    # Sumatra Utara
    geom_sf(data = indo_sumut,
            fill = "#cde8b6",   
            color = "grey70",
            linewidth = 0.2) +
    
    # Heatmap
    geom_tile(data = data_plot,
              aes(x = x, y = y, fill = var1.pred),
              width = 0.05,
              height = 0.05,
              alpha = 0.75) +
    
    # Wilayah penelitian
    geom_sf(data = region_select,
            fill = NA,
            color = "black",
            linewidth = 0.4) +
    
    # Stasiun point
    geom_point(data = data_spasial,
               aes(x = Longitude, y = Latitude),
               color = "red",
               size = 2) +
    geom_text(data = data_spasial,
              aes(x = Longitude, y = Latitude, label = Station),
              size = 3,
              hjust = -0.1,
              vjust = -0.5) +
    
    scale_fill_gradientn(
      name = "Return Level (mm)",
      colours = c("#0d0887", "#6a00a8", "#b12a90", "#e16462", "#fca636", "#f0f921"),
      values = scales::rescale(c(140, 180, 220, 260, 300, 340)),
      limits = c(140, 340),
      breaks = seq(140, 340, by = 40),
      oob = scales::squish
    ) +
    
    coord_sf(expand = FALSE) +
    
    labs(
      title = title,
      x = "Longitude",
      y = "Latitude"
    ) +
    
    theme_minimal() +
    theme(
      # Sea
      panel.background = element_rect(fill = "#b9e3f9", color = NA),
      # Grid
      panel.grid.major = element_line(color = "white", linewidth = 0.1),
      # Title
      plot.title = element_text(hjust = 0.5)
    )
}

plot_10 <- plot_rl(idwr_10, "Return Level 10 Year (Fused Ridge)")
plot_10
plot_25 <- plot_rl(idwr_25, "Return Level 25 Year (Fused Ridge)")
plot_25
plot_100 <- plot_rl(idwr_100, "Return Level 100 Year (Fused Ridge)")
plot_100
