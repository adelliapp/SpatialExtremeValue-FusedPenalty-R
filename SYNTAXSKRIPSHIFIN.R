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
library(evir)
library(ggspatial)
library(sf)
library(akima)       # interpolasi IDW
library(splines)
library(viridis)
library(goftest)
library(geosphere)   # distm/jarak euclidian
library(GGally)      # plot korelasi antar lokasi
library(nortest)     # uji A-D
library(geostatsp)   # F-Madogram
library(reshape2)    # Gaussian Copula
library(corrplot)    # correlation location 
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



# =========================
#      DESKRIPSI DATA
# =========================
# Statistik deskriptif
summary_data <- data %>%
  group_by(Stasiun) %>%  # Kelompokkan berdasarkan Lokasi
  summarise(
    Min = min(Curah.Hujan, na.rm = TRUE),
    Q1 = quantile(Curah.Hujan, 0.25, na.rm = TRUE),
    Median = median(Curah.Hujan, na.rm = TRUE),
    Mean = mean(Curah.Hujan, na.rm = TRUE),
    Q3 = quantile(Curah.Hujan, 0.75, na.rm = TRUE),
    Max = max(Curah.Hujan, na.rm = TRUE),
    SD  = sd(Curah.Hujan, na.rm = TRUE)
  )
print(summary_data)

# Boxplot Sebaran Curah Hujan di Tiap Stasiun
ggplot(data, aes(x = Stasiun, y = Curah.Hujan, fill = Stasiun)) +
  geom_boxplot(alpha = 0.7, outlier.color = "red") +
  labs(
    title = "Sebaran Curah Hujan per Stasiun",
    x = "Stasiun Pengamatan",
    y = "Curah Hujan (mm)"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# Titik Sebaran Curah Hujan di Tiap Stasiun
ggplot(data, aes(x = Stasiun, y = Curah.Hujan, color = Stasiun)) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +  # titik sebar dengan sedikit noise horizontal
  labs(
    title = "Sebaran Titik Curah Hujan per Stasiun",
    x = "Stasiun Pengamatan",
    y = "Curah Hujan (mm)"
  ) +
  theme_minimal() +
  theme(legend.position = "none")





# =========================
#   KORELASI ANTAR LOKASI
# =========================
# Ubah ke format wide (tiap stasiun jadi kolom)
data_wide <- read.csv("D:/DontTouch!/KULEEYEEAHH SMT 7/CMEMEW/PENTING/Korelasi Data Curah Hujan Sumut 14-24.csv")
summary(data_wide)

# Ambil hanya kolom stasiun (hilangkan 'id')
data_plot <- data_wide %>% dplyr::select(-id)
# Hapus baris dengan NA 
data_plot <- na.omit(data_plot)

# Jumlah observasi yang bisa dianalisis
nrow(data_plot)

# Plot korelasi antar stasiun
ggpairs(data_plot,
        title = "Plot Korelasi Curah Hujan Ekstrem Antar Stasiun",
        upper = list(continuous = GGally::ggally_cor),
        lower = list(continuous = GGally::ggally_points),
        diag = list(continuous = GGally::ggally_densityDiag))

# Korelasi dalam bentuk tabel
cor(data_plot, use = "complete.obs")




# ====================================
#   KORELASI ANTAR DATA TIAP STASIUN
# ====================================
# Fungsi Ljung-Box
uji_dependensi_stasiun <- function(data, nama_stasiun, lag_val = 5) {
  x <- data$Curah.Hujan  # data asli (time series)
  hasil <- list(
    Stasiun = nama_stasiun,
    n_data = length(x)
  )
  
  if (length(x) > lag_val) {
    lb_test <- Box.test(x, lag = lag_val, type = "Ljung-Box")
    hasil$lag <- lag_val
    hasil$ljung_pval <- lb_test$p.value
    
    # Kesimpulan
    if (lb_test$p.value > 0.05) {
      hasil$kesimpulan <- "Independen"
    } else {
      hasil$kesimpulan <- "Ada ketergantungan"
    }
    
  } else {
    hasil$lag <- lag_val
    hasil$ljung_pval <- NA
    hasil$kesimpulan <- "Data tidak cukup"
  }
  
  return(hasil)
}

# Loop semua stasiun
hasil_dependensi <- list()
stasiun_list <- unique(data$Stasiun)
for (s in stasiun_list) {
  df_s <- data %>% dplyr::filter(Stasiun == s)
  hasil <- uji_dependensi_stasiun(df_s, s, lag_val = 5)
  hasil_dependensi[[s]] <- hasil
}

# Gabungkan hasil jadi data.frame
hasil_dependensi_df <- do.call(rbind, lapply(hasil_dependensi, as.data.frame))
print(hasil_dependensi_df)

# Visualisasi plot
ggplot(hasil_dependensi_df, aes(x = Stasiun, y = ljung_pval)) +
  geom_point(size = 4, color = "blue") +
  geom_hline(yintercept = 0.05, color = "red", linetype = "dashed") +
  labs(title = "Uji Ljung-Box per Stasiun",
       y = "p-value", x = "Stasiun") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




# ========================
#   PEMILIHAN THRESHOLD
# ========================
# Metode Dinamis Threshold 
estimasi_gpd_per_lokasi <- function(df, lokasi) {
  data_hujan <- df %>% filter(Stasiun == lokasi) %>% pull(Curah.Hujan)
  percentiles <- seq(0.50, 0.95, by = 0.01)
  thresholds <- quantile(data_hujan, percentiles, na.rm = TRUE)
  
  results <- data.frame()
  
  for (i in seq_along(thresholds)) {
    u <- thresholds[i]
    if (sum(data_hujan > u) < 10) next
    
    fit_gpd <- ismev::gpd.fit(data_hujan, threshold = u, show = FALSE)
    if (is.null(fit_gpd)) next
    
    data_ekstrem <- fit_gpd$data[fit_gpd$data > u] - u
    
    xi <- fit_gpd$mle[2]
    sigma <- fit_gpd$mle[1]
    data_norm <- data_ekstrem / sigma
    
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
    
    n_exc <- length(data_ekstrem)
    
    results <- rbind(results, data.frame(
      Percentile = percentiles[i] * 100,
      Threshold = round(u, 2),
      xi = round(xi, 2),
      sigma = round(sigma, 2),
      AD_Stat = round(ad$statistic, 3),
      AD_p_value = round(ad$p.value, 3),
      Jumlah_Exceedance = n_exc,
      Stasiun = lokasi
    ))
  }
  
  results <- results %>% arrange(Threshold)
  results <- results %>% mutate(
    MRL = map_dbl(Threshold, function(u) {
      exc <- data_hujan[data_hujan > u]
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
    mutate(Threshold_Terbaik = TRUE)
  
  results <- results %>%
    mutate(Threshold_Terbaik = Threshold %in% best$Threshold)
  
  return(results)
}

# MRL Plot untuk semua Stasiun
mrl_all <- map_df(unique(data$Stasiun), function(stasiun_nama) {
  data_stasiun <- data %>% filter(Stasiun == stasiun_nama)
  
  th_seq <- quantile(data_stasiun$Curah.Hujan,
                     probs = seq(0.50, 0.95, 0.01),
                     na.rm = TRUE)
  
  mrl <- map_dbl(th_seq, function(u) {
    exc <- data_stasiun$Curah.Hujan[data_stasiun$Curah.Hujan > u]
    if (length(exc) > 0) mean(exc - u) else NA_real_
  })
  
  data.frame(
    Stasiun = stasiun_nama,
    Threshold = th_seq,
    MRL = mrl
  )
})

threshold_summary <- unique(data$Stasiun) %>%
  map_df(~estimasi_gpd_per_lokasi(data, .x)) %>%
  dplyr::filter(Threshold_Terbaik)

mrl_all <- mrl_all %>%
  left_join(
    threshold_summary %>% dplyr::select(Stasiun, Threshold),
    by = "Stasiun",
    suffix = c("", "_final")
  )

ggplot(mrl_all, aes(x = Threshold, y = MRL)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(color = "steelblue", size = 1.5) +
  
  # garis threshold terbaik
  geom_vline(aes(xintercept = Threshold_final),
             linetype = "dashed",
             color = "red",
             linewidth = 0.8) +
  facet_wrap(~Stasiun, scales = "free") +
  labs(
    title = "Mean Residual Life (MRL) Plot untuk Seluruh Stasiun",
    x = "Threshold (mm)",
    y = "Mean Excess"
  ) +
  theme_minimal(base_size = 13)

# Hasil Threshold Tiap Stasiun
print(threshold_summary)
write.csv(threshold_summary, "D:/DontTouch!/KULEEYEEAHH SMT 7/CMEMEW/PENTING/CRHJN_Threshold.csv", row.names = FALSE)





# ====================
#   DATA EXCEEDANCE
# ====================
data_exceedance <- data %>%
  left_join(threshold_summary %>% 
              rename(Threshold = Threshold) %>% 
              dplyr::select(Stasiun, Threshold), 
            by = "Stasiun") %>%
  filter(Curah.Hujan > Threshold) %>%  # hanya data ekstrem
  mutate(
    Excess = Curah.Hujan - Threshold
  ) %>%
  dplyr::select(Stasiun, Curah.Hujan, Threshold, Excess, Latitude, Longitude)
head(data_exceedance)

data_exceedance %>%
  group_by(Stasiun) %>%
  summarise(n_exceed = n(), .groups = "drop")

hist(data_exceedance$Excess, breaks = 30)
write.csv(data_exceedance, "D:/DontTouch!/KULEEYEEAHH SMT 7/CMEMEW/PENTING/CRHJN_data_exceedances.csv", row.names = FALSE)





# =============================
#   FUNGSI LOG-LIKELIHOOD GPD 
# =============================
# Parameter bentuk = konstan
xi_data <- data_exceedance$Excess
fit_xi <- ismev::gpd.fit(xi_data, threshold = 0, show = FALSE)
xi_konstan <- fit_xi$mle[2]
print(xi_konstan)

data_list <- split(data_exceedance$Excess,
                   data_exceedance$Stasiun)
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
# HITUNG MATRIKS BOBOT (Wpq)
# Ambil koordinat unik per lokasi
lokasi_coords <- data_exceedance %>%
  group_by(Stasiun) %>%
  summarise(
    Latitude = mean(Latitude, na.rm = TRUE),
    Longitude = mean(Longitude, na.rm = TRUE)
  )

# Hitung jarak Haversine antar lokasi (meter)
coords_matrix <- as.matrix(lokasi_coords[, c("Longitude", "Latitude")])
jarak_matrix <- distm(coords_matrix, fun = distHaversine)

library(igraph)
g <- graph_from_adjacency_matrix(
  jarak_matrix,
  mode = "undirected",
  weighted = TRUE,
  diag = FALSE
)
V(g)$name <- lokasi_coords$Stasiun              # beri nama vertex
mst_g <- mst(g, weights = E(g)$weight)          # cari MST
edge_df <- as_data_frame(mst_g, what = "edges") # ambil pasangan tetangga

# Ubah ke bobot invers (Wpq = 1 / jarak)
Wpq_matrix <- 1 / jarak_matrix
diag(Wpq_matrix) <- 0  # nolkan diagonal agar tidak self-loop
Wpq_matrix_scaled <- Wpq_matrix / mean(Wpq_matrix[Wpq_matrix > 0])
edge_df$Wpq <- mapply(
  function(a, b){
    i <- match(a, lokasi_coords$Stasiun)
    j <- match(b, lokasi_coords$Stasiun)
    
    Wpq_matrix_scaled[i, j]
  },
  edge_df$from,
  edge_df$to
)
edge_df <- edge_df %>%
  rename(
    Lokasi_1 = from,
    Lokasi_2 = to
  )
print(edge_df)

# Simpan hasil ke Struktur Siap Input
# Buat list data per lokasi (exceedance)
data_list <- data_exceedance %>%
  group_by(Stasiun) %>%
  summarise(Excess = list(Excess), .groups = "drop")

# Susun jadi list input
input_fused_penalty <- list(
  xi = xi_konstan,
  lokasi = lokasi_coords,
  edge = edge_df,
  W_matrix = Wpq_matrix_scaled,
  excess_per_lokasi = data_list
)

saveRDS(input_fused_penalty, "input_fused_penalty.rds")  # jika ingin disimpan ke file



# LOAD DATA INPUT (berisi: xi, lokasi, edge, W_matrix, excess_per_lokasi)
input <- readRDS("input_fused_penalty.rds")

# Ambil masing-masing komponen
data_list <- input$excess_per_lokasi$Excess   # list excess per stasiun
n_lokasi <- length(data_list)
stasiun_names <- input$lokasi$Stasiun
xi_konstan <- input$xi
edge_df <- input$edge

# Inisialisasi Sigma pakai estimasi sigma tanpa penalty
estimasi_gpd_tanpa_regularisasi <- function(df, lokasi){
  data_excess <- df %>%
    filter(Stasiun == lokasi)
  
  excess <- data_excess$Excess
  
  fit <- ismev::gpd.fit(
    excess,
    threshold = 0,
    show = FALSE
  )
  
  data.frame(
    Stasiun = lokasi,
    Sigma = fit$mle[1],
    Xi = fit$mle[2]
  )
}
stasiun_list <- unique(data_exceedance$Stasiun)
hasil_estimasi <- map_dfr(stasiun_list, ~ estimasi_gpd_tanpa_regularisasi(data_exceedance, .x))
print(hasil_estimasi)
sigma_init <- input$excess_per_lokasi %>%
  left_join(hasil_estimasi, by = "Stasiun") %>%
  pull(Sigma)
identical(input$excess_per_lokasi$Stasiun, 
          hasil_estimasi$Stasiun[match(input$excess_per_lokasi$Stasiun, hasil_estimasi$Stasiun)])   # Harus TRUE

print(input_fused_penalty)



# ===============
#   FUSED LASSO
# ===============
# Penalti Fused Lasso
penalty_lasso <- function(sigma, edge_df, stasiun_names, lambda){
  penalty <- 0
  
  for(k in 1:nrow(edge_df)){
    p <- match(edge_df$Lokasi_1[k], stasiun_names)
    q <- match(edge_df$Lokasi_2[k], stasiun_names)
    
    penalty <- penalty +
      edge_df$Wpq[k] *
      abs(sigma[p] - sigma[q])
  }
  
  return(lambda * penalty)
}

# Fungsi objective lasso untuk optim
objective_lasso <- function(sigma,
                            data_list,
                            xi,
                            edge_df,
                            stasiun_names,
                            lambda){
  
  ll <- loglik_gpd(sigma, xi, data_list)
  
  pen_lasso <- penalty_lasso(
    sigma,
    edge_df,
    stasiun_names,
    lambda
  )
  
  return(-(ll - pen_lasso))
}

# Seleksi Lamda (Parameter Regularisasi)
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
    lower = rep(1e-6, n_lokasi),
    data_list = data_list,
    xi = xi_konstan,
    edge_df = edge_df,
    stasiun_names = stasiun_names,
    lambda = lambda
  )
  
  sigma_tmp <- res$par
  loglik_tmp <- -objective_lasso(
    sigma_tmp,
    data_list,
    xi_konstan,
    edge_df,
    stasiun_names,
    lambda
  )
  hess_tmp <- hessian(
    func = objective_lasso,
    x = sigma_tmp,
    data_list = data_list,
    xi = xi_konstan,
    edge_df = edge_df,
    stasiun_names = stasiun_names,
    lambda = lambda
  )
  grad_tmp <- grad(
    func = objective_lasso,
    x = sigma_tmp,
    data_list = data_list,
    xi = xi_konstan,
    edge_df = edge_df,
    stasiun_names = stasiun_names,
    lambda = lambda
  )
  H_inv <- tryCatch(solve(hess_tmp), error = function(e) MASS::ginv(hess_tmp))
  J_matrix <- grad_tmp %*% t(grad_tmp)
  
  tic_lambda_lasso$TIC_lambda[i] <- -2 * (loglik_tmp - sum(diag(H_inv %*% J_matrix)))
}

best_lambda_lasso <- tic_lambda_lasso$lambda[which.min(tic_lambda_lasso$TIC_lambda)]
print(best_lambda_lasso)

# Estimasi final dengan lambda terbaik
res_best_lasso <- optim(
  par = sigma_init,
  fn = objective_lasso,
  method = "L-BFGS-B",
  lower = rep(1e-6, n_lokasi),
  data_list = data_list,
  xi = xi_konstan,
  edge_df = edge_df,
  stasiun_names = stasiun_names,
  lambda = best_lambda_lasso,
  control = list(maxit = 5000, factr = 1e-10)
)

# Hasil Fused Lasso
sigma_hat_lasso <- res_best_lasso$par
hasil_lasso <- data.frame(
  Stasiun = stasiun_names,
  Sigma_FusedLasso = sigma_hat_lasso
)
print(hasil_lasso)



# ===============
#   FUSED RIDGE
# ===============
# Penalti Fused Ridge
penalty_ridge <- function(sigma, edge_df, stasiun_names, lambda){
  penalty <- 0
  
  for(k in 1:nrow(edge_df)){
    p <- match(edge_df$Lokasi_1[k], stasiun_names)
    q <- match(edge_df$Lokasi_2[k], stasiun_names)
    
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
                            stasiun_names,
                            lambda){
  
  ll <- loglik_gpd(sigma, xi, data_list)
  
  pen_ridge <- penalty_ridge(
    sigma,
    edge_df,
    stasiun_names,
    lambda
  )
  
  return(-(ll - pen_ridge))
}

# Seleksi Lamda (Parameter Regularisasi)
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
    stasiun_names = stasiun_names,
    lambda = lambda
  )
  
  sigma_tmp <- res$par
  loglik_tmp <- -objective_ridge(
    sigma_tmp,
    data_list,
    xi_konstan,
    edge_df,
    stasiun_names,
    lambda
  )
  hess_tmp <- hessian(
    func = objective_ridge,
    x = sigma_tmp,
    data_list = data_list,
    xi = xi_konstan,
    edge_df = edge_df,
    stasiun_names = stasiun_names,
    lambda = lambda
  )
  grad_tmp <- grad(
    func = objective_ridge,
    x = sigma_tmp,
    data_list = data_list,
    xi = xi_konstan,
    edge_df = edge_df,
    stasiun_names = stasiun_names,
    lambda = lambda
  )
  H_inv <- tryCatch(solve(hess_tmp), error = function(e) MASS::ginv(hess_tmp))
  J_matrix <- grad_tmp %*% t(grad_tmp)
  
  tic_lambda_ridge$TIC_lambda[i] <- -2 * (loglik_tmp - sum(diag(H_inv %*% J_matrix)))
}

best_lambda_ridge <- tic_lambda_ridge$lambda[which.min(tic_lambda_ridge$TIC_lambda)]
print(best_lambda_ridge)

# Estimasi final dengan lambda terbaik
res_best_ridge <- optim(
  par = sigma_init,
  fn = objective_ridge,
  method = "L-BFGS-B",
  lower = rep(1e-6, n_lokasi),
  data_list = data_list,
  xi = xi_konstan,
  edge_df = edge_df,
  stasiun_names = stasiun_names,
  lambda = best_lambda_ridge,
  control = list(maxit = 5000, factr = 1e-10)
)

# Hasil Fused Ridge
sigma_hat_ridge <- res_best_ridge$par
hasil_ridge <- data.frame(
  Stasiun = stasiun_names,
  Sigma_FusedRidge = sigma_hat_ridge
)
print(hasil_ridge)





# ===========================
#   PLOT HASIL REGULARISASI
# ===========================
plot(
  hasil_lasso$Sigma_FusedLasso,
  hasil_ridge$Sigma_FusedRidge,
  xlab = "Fused Lasso",
  ylab = "Fused Ridge",
  main = "Perbandingan Sigma",
  pch = 19
)

abline(a = 0, b = 1, col = "red", lwd = 2)




# ===========
#   QQ PLOT 
# ===========
# Fused Lasso
qgpd_theoretical <- function(p, xi, beta) {
  if (xi == 0) return(-beta * log(1 - p))
  return(beta * ((1 - p)^(-xi) - 1) / xi)
}
par(mfrow = c(1, length(data_list)), oma = c(0, 0, 3, 0))  # ruang atas untuk judul
for (i in seq_along(data_list)) {
  excess_data <- sort(data_list[[i]])
  n <- length(excess_data)
  p_empirical <- ppoints(n)
  
  beta_i <- hasil_lasso$Sigma_FusedLasso[i]
  
  q_theoretical <- qgpd_theoretical(p_empirical, xi_konstan, beta_i)
  
  plot(q_theoretical, excess_data,
       main = hasil_ridge$Stasiun[i], 
       cex.main = 1.2,
       xlab = "Theoretical Quantiles",
       ylab = "Empirical Quantiles")
  
  abline(0, 1, col = "red", lwd = 2)
}
title(main = "QQ Plot GPD dengan Estimasi MLE (Regularisasi Fused Lasso)",
      outer = TRUE, cex.main = 1.5)

# Fused Ridge
qgpd_theoretical <- function(p, xi, beta) {
  if (xi == 0) return(-beta * log(1 - p))
  return(beta * ((1 - p)^(-xi) - 1) / xi)
}
par(mfrow = c(1, length(data_list)), oma = c(0, 0, 3, 0))    # ruang atas untuk judul
for (i in seq_along(data_list)) {
  excess_data <- sort(data_list[[i]])
  n <- length(excess_data)
  p_empirical <- ppoints(n)
  
  beta_i <- hasil_ridge$Sigma_FusedRidge[i]
  
  q_theoretical <- qgpd_theoretical(p_empirical, xi_konstan, beta_i)
  
  plot(q_theoretical, excess_data,
       main = hasil_ridge$Stasiun[i], 
       cex.main = 1.2,
       xlab = "Theoretical Quantiles",
       ylab = "Empirical Quantiles")
  
  abline(0, 1, col = "red", lwd = 2)
}
title(main = "QQ Plot GPD dengan Estimasi MLE (Regularisasi Fused Ridge)",
      outer = TRUE, cex.main = 1.5)




# ==============================
#   UJI KESESUAIAN DISTRIBUSI
# ==============================
ad_gpd_custom <- function(data_excess, xi, beta) {
  n <- length(data_excess)
  # Transformasi ke Uniform (0,1)
  u <- evir::pgpd(data_excess, xi = xi, beta = beta)
  u <- sort(u[is.finite(u) & u > 0 & u < 1])  # Pastikan dalam (0,1)
  n <- length(u)
  
  i <- 1:n
  A2 <- -n - (1/n) * sum((2*i - 1) * (log(u) + log(1 - rev(u))))
  
  return(list(stat = A2))
}

uji_kesesuaian_gpd <- function(data_ekstrem, sigma, xi_hat) {
  # Anderson-Darling Khusus GPD
  ad_result <- ad_gpd_custom(data_ekstrem, xi = xi_hat, beta = sigma)
  AD_Stat <- ad_result$stat
  AD_Kesimpulan <- ad_result$kesimpulan
  
  # Likelihood Ratio Test (uji eksponensial vs GPD)
  n <- length(data_ekstrem)
  dens <- evir::dgpd(data_ekstrem, xi = xi_hat, beta = sigma)
  dens <- dens[is.finite(dens) & dens > 0]
  loglik_gpd <- sum(log(dens))
  scale_exp <- mean(data_ekstrem)
  loglik_exp <- sum(dexp(data_ekstrem, rate = 1 / scale_exp, log = TRUE))
  LRT_Stat <- -2 * (loglik_exp - loglik_gpd)
  Chi_Crit <- qchisq(0.95, df = 1)
  LRT_Kesimpulan <- ifelse(LRT_Stat > Chi_Crit, "Tolak H0", "Gagal Tolak H0")
  
  return(data.frame(
    AD_Stat = round(AD_Stat, 4),
    LRT_Stat = round(LRT_Stat, 4),
    Chi_Crit = round(Chi_Crit, 4),
    LRT_Kesimpulan = LRT_Kesimpulan
  ))
}

# Susun ulang nama2 untuk sinkronisasi
stasiun_list <- unique(data_exceedance$Stasiun)
names(data_list) <- input$excess_per_lokasi$Stasiun  # untuk hasil fused

# Tambahkan hasil uji kesesuaian untuk setiap hasil regularisasi
# Fungsi bantu gabungan uji per lokasi dan jenis regularisasi
uji_per_lokasi <- function(stasiun, sigma_vec, label) {
  idx <- which(names(data_list) == stasiun)
  excess_data <- data_list[[idx]]
  sigma_val <- sigma_vec[idx]
  hasil_uji <- uji_kesesuaian_gpd(excess_data, sigma_val, xi_konstan)
  hasil_uji$Stasiun <- stasiun  # <--- tambahkan ini
  colnames(hasil_uji)[colnames(hasil_uji) != "Stasiun"] <- 
    paste0(colnames(hasil_uji)[colnames(hasil_uji) != "Stasiun"], "_", label)
  return(hasil_uji)
}

# Siapkan hasil uji AD dan LRT untuk Fused Lasso
hasil_uji_lasso <- map2_dfr(
  data_list, hasil_lasso$Sigma_FusedLasso,
  ~ uji_kesesuaian_gpd(.x, .y, xi_konstan)
) %>%
  mutate(Stasiun = hasil_lasso$Stasiun,
         Metode = "Fused Lasso") %>%
  dplyr::select(Stasiun, Metode, everything())

# Siapkan hasil uji AD dan LRT untuk Fused Ridge
hasil_uji_ridge <- map2_dfr(
  data_list, hasil_ridge$Sigma_FusedRidge,
  ~ uji_kesesuaian_gpd(.x, .y, xi_konstan)
) %>%
  mutate(Stasiun = hasil_ridge$Stasiun,
         Metode = "Fused Ridge") %>%
  dplyr::select(Stasiun, Metode, everything())

# Gabungkan hasil kedua metode
hasil_uji_gpd <- bind_rows(hasil_uji_lasso, hasil_uji_ridge)

# Simpan Hasil Estimasi dan Uji Kesesuaian Distribusi
hasil_uji_lasso$Sigma_Hat <- sigma_hat_lasso
hasil_uji_ridge$Sigma_Hat <- sigma_hat_ridge
hasil_uji_gpd <- bind_rows(hasil_uji_lasso, hasil_uji_ridge)
hasil_uji_gpd$Xi_Konstan <- xi_konstan
print(hasil_uji_gpd)
write.csv(hasil_uji_gpd, "D:/DontTouch!/KULEEYEEAHH SMT 7/CMEMEW/PENTING/CRHJN_hasil_estimasi_dan_uji_distribusi.csv", row.names = FALSE)





# =============================
#   UJI DEPENDENSI F MADOGRAM
# =============================
data_wide <- read.csv("D:/DontTouch!/KULEEYEEAHH SMT 7/CMEMEW/PENTING/Korelasi Data Curah Hujan Sumut 14-24.csv")
summary(data_wide)

# Ambil data Latitude dan Longitude
lokasi <- data %>%
  dplyr::select(Stasiun, Longitude, Latitude) %>%
  dplyr::distinct() %>%          # Ambil satu baris unik per stasiun (karena koordinat sama terus)
  dplyr::mutate(Stasiun = gsub(" ", ".", Stasiun))    # samakan format nama

# Fungsi Hitung Koefisien Ekstermal pakai F-Madogram
hitung_fmadogram <- function(data_wide, lokasi) {
  stasiun <- setdiff(names(data_wide), "id")       # ambil nama stasiun dari kolom
  
  lokasi_pair <- combn(stasiun, 2, simplify = FALSE)
  
  map_dfr(lokasi_pair, function(pair) {
    loc1 <- pair[1]; loc2 <- pair[2]
    x <- data_wide[[loc1]]
    y <- data_wide[[loc2]]
    
    # hapus NA
    ok <- complete.cases(x, y)
    x <- x[ok]; y <- y[ok]
    n <- length(x)
    
    if (n < 30) return(NULL)
    
    # pseudo-uniform pakai rank
    u1 <- rank(x, ties.method="average") / (n+1)
    u2 <- rank(y, ties.method="average") / (n+1)
    
    # F-madogram
    nu_hat <- 0.5 * mean(abs(u1 - u2))
    
    # Extremal coefficient
    theta_hat <- ifelse(nu_hat >= 1/6, 2,
                        (1 + 2*nu_hat) / (1 - 2*nu_hat))
    
    # koordinat (FIX urutan!)
    coords <- lokasi %>%
      filter(Stasiun %in% pair) %>%
      arrange(match(Stasiun, pair)) %>%
      dplyr::select(Longitude, Latitude)
    if(nrow(coords) < 2) return(NULL)
    d <- geosphere::distHaversine(as.matrix(coords[1, ]), as.matrix(coords[2, ]))
    
    tibble(
      Lokasi1 = loc1,
      Lokasi2 = loc2,
      Distance_km = d / 1000,
      Fmadogram = nu_hat,
      Extremal_Coeff = theta_hat
    )
  })
}

fmadogram_result <- hitung_fmadogram(data_wide, lokasi)
print(fmadogram_result)

# Plot F Madogram
par(mfrow = c(1,1))   # balikin ke 1 plot
plot(fmadogram_result$Distance_km, 
     fmadogram_result$Extremal_Coeff,
     pch = 19,
     col = "blue",
     xlab = "Jarak (km)",
     ylab = "Extremal Coefficient",
     main = "F-Madogram Plot")





# ==========================================
#   ESTIMASI PARAMETER GAUSSIAN COPULA (R)
# ==========================================
fix_station <- function(x) {
  x %>%
    str_replace_all(" ", "\\.") %>%  # spasi -> titik
    stringr::str_squish()            # rapikan spasi
}
threshold_summary <- threshold_summary %>%
  dplyr::mutate(Stasiun = fix_station(Stasiun))
hasil_uji_gpd <- hasil_uji_gpd %>%
  dplyr::mutate(Stasiun = fix_station(Stasiun))
data_exceedance <- data_exceedance %>%
  dplyr::mutate(Stasiun = fix_station(Stasiun))

transform_uniform_semiparametric <- function(data_wide,
                                             parameter_gpd,
                                             threshold_summary){
  hasil <- data_wide
  nama_stasiun <- colnames(data_wide)[-1]
  
  for(st in nama_stasiun){
    x <- data_wide[[st]]
    threshold <- threshold_summary %>%
      filter(Stasiun == st) %>%
      pull(Threshold)
    sigma <- parameter_gpd %>%
      filter(Stasiun == st) %>%
      pull(Sigma_Hat)
    xi <- parameter_gpd %>%
      filter(Stasiun == st) %>%
      pull(Xi_Konstan)
    
    # proporsi exceedance
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
    
    hasil[[st]] <- u
  }
  
  return(hasil)
}

parameter_lasso <- hasil_uji_gpd %>%
  filter(Metode == "Fused Lasso")
data_uniform_lasso <- transform_uniform_semiparametric(
  data_wide,
  parameter_lasso,
  threshold_summary)

parameter_ridge <- hasil_uji_gpd %>%
  filter(Metode == "Fused Ridge")
data_uniform_ridge <- transform_uniform_semiparametric(
  data_wide,
  parameter_ridge,
  threshold_summary)

# Fungsi Gaussian Copula
gauss_copula <- function(data_uniform, nama_metode){
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

# Simpan hasil copula Gaussian
saveRDS(copula_lasso, "copula_lasso.rds")
saveRDS(copula_ridge, "copula_ridge.rds")





# ================
#   RETURN LEVEL
# ================
# Fungsi Return Level Uni
hitung_return_level_uni <- function(exceedance, 
                                    parameter_gpd,
                                    threshold,
                                    T_tahun){
  
  T_hari <- T_tahun * 365
  
  lambda_table <- exceedance %>%
    group_by(Stasiun) %>%
    summarise(lambda_u = n()/4018, .groups = "drop")
  
  param <- parameter_gpd %>%
    dplyr::select(Stasiun, Sigma_Hat, Xi_Konstan) %>%
    left_join(threshold %>%
                dplyr::select(Stasiun, Threshold),
              by = "Stasiun") %>%
    left_join(lambda_table, by="Stasiun") 
  
  probabilitas <- 1 - 1/(T_hari * param$lambda_u)
  
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

#  Fungsi Invers Semi-Parametric
inverse_semiparametric <- function(u, x_data,
                                   threshold,
                                   sigma,
                                   xi){
  
  lambda_u <- mean(x_data > threshold)
  if(u <= (1-lambda_u)){
    x_bawah <- x_data[x_data <= threshold]
    return(
      quantile(
        x_bawah,
        probs = u/(1-lambda_u),
        type = 8))
  }else{
    p_tail <- (u-(1-lambda_u))/lambda_u
    if(abs(xi) < 1e-6){
      return(
        threshold - sigma*log(1-p_tail))
    }else{
      return(
        threshold +
          (sigma/xi)*
          ((1-p_tail)^(-xi)-1))
    }
  }
}

# Fungsi Return Level Mix (dengan periode dalam bentuk hari)
hitung_return_level_multi <- function(data_wide, 
                                      parameter_gpd,
                                      threshold,
                                      R_hat,
                                      copula_model,
                                      T_tahun, 
                                      n_simulasi = 50000) {
  param <- parameter_gpd %>%
    dplyr::select(Stasiun, Sigma_Hat, Xi_Konstan) %>%
    left_join(
      threshold %>%
        dplyr::select(Stasiun, Threshold),
      by = "Stasiun") %>%
    arrange(match(Stasiun, colnames(R_hat)))
  
  d <- ncol(R_hat)
  
  # Simulasi Gaussian Copula
  x_sim <- matrix(NA, nrow = n_simulasi, ncol = d)
  u_sim <- rCopula(n_simulasi, copula_model)
  
  for(i in 1:d){
    stasiun <- param$Stasiun[i]
    sigma   <- param$Sigma_Hat[i]
    xi      <- param$Xi_Konstan[i]
    u_thr   <- param$Threshold[i]
    
    x_data <- data_wide[[stasiun]]
    x_sim[,i] <- sapply(
      u_sim[,i],
      inverse_semiparametric,
      x_data = x_data,
      threshold = u_thr,
      sigma = sigma,
      xi = xi)
  }
  
  system_extreme <- apply(x_sim, 1, max)
  
  return_multi <- quantile(system_extreme, probs = 1 - 1/(365*T_tahun))
  
  return(as.numeric(return_multi))
}

# Periode return level
periode_list <- c(10, 25, 100)

# List kosong untuk simpan hasil
hasil_RL <- list()

set.seed(123)
for (T in periode_list) {
  RL_Uni_Lasso <- hitung_return_level_uni(data_exceedance, 
                                          parameter_lasso,
                                          threshold_summary,
                                          T_tahun = T)
  RL_Uni_Ridge <- hitung_return_level_uni(data_exceedance,
                                          parameter_ridge,
                                          threshold_summary,
                                          T_tahun = T)
  
  RL_Multi_Lasso <- hitung_return_level_multi(data_wide,
                                              parameter_lasso,
                                              threshold_summary,
                                              R_hat_lasso,
                                              cop_model_lasso,
                                              T_tahun = T)
  
  RL_Multi_Ridge <- hitung_return_level_multi(data_wide,
                                              parameter_ridge,
                                              threshold_summary,
                                              R_hat_ridge,
                                              cop_model_ridge,
                                              T_tahun = T)
  
  # Gabungkan univariat
  df_temp <- RL_Uni_Lasso %>%
    dplyr::select(Stasiun, RL_Lasso = Return_Level) %>%
    left_join(
      RL_Uni_Ridge %>%
        dplyr::select(Stasiun, RL_Ridge = Return_Level),
      by = "Stasiun"
    )
  
  # Tambahkan multivariat (sama untuk semua stasiun → jadi kolom)
  df_temp$RL_System_Lasso <- RL_Multi_Lasso
  df_temp$RL_System_Ridge <- RL_Multi_Ridge
  
  # Tambahkan info periode
  df_temp$Periode_Tahun <- T
  
  # Simpan ke list
  hasil_RL[[as.character(T)]] <- df_temp
}

# Gabungkan semua periode
return_level_final <- bind_rows(hasil_RL)

# Lihat hasil
print(return_level_final)

write.csv(return_level_final, "D:/DontTouch!/KULEEYEEAHH SMT 7/CMEMEW/PENTING/CRHJN_return_level.csv", row.names = FALSE)





# ===================================
#   PERBANDINGAN REGULARISASI (TIC)
# ===================================
# Fungsi LogLik (PDF) GPD total
loglik_gpd <- function(x, sigma, xi) {
  # Cek kelayakan nilai sigma
  if (any(sigma <= 0)) return(-Inf)
  
  # Cek syarat batas GPD: 1 + xi * x / sigma harus > 0
  cond <- 1 + xi * x / sigma
  if (any(cond <= 0)) return(-Inf)
  if (abs(xi) < 1e-6) {
    ll <- -sum(log(sigma) + x / sigma)
  } else {
    ll <- sum(-log(sigma) - (1/xi + 1) * log(cond))
  }
  return(ll)
}

# Log-likelihood per observasi
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

# Fungsi Hitung TIC Untuk Model Uni
hitung_TIC_Uni <- function(data, sigma_vector, xi) {
  stasiun_unik <- unique(data$Stasiun)
  n_stasiun <- length(sigma_vector)
  
  # Mapping stasiun ke indeks angka (1 sampai 5) untuk presisi looping
  data$Stasiun_index <- match(data$Stasiun, stasiun_unik)
  
  # sigma untuk tiap observasi
  sigma_obs <- sigma_vector[data$Stasiun_index]
  
  # log-likelihood total
  loglik <- loglik_gpd(
    x = data$Excess,
    sigma = sigma_obs,
    xi = xi
  )
  
  # Gabungkan parameter (5 sigma + 1 xi)
  param <- c(sigma_vector, xi)
  
  # Hitung Hessian (J) menggunakan Likelihood Total
  grad_func <- function(p){
    sigma_temp <- p[1:n_stasiun]
    xi_temp <- p[length(p)]
    sigma_obs_temp <- sigma_temp[data$Stasiun_index]
    
    return(loglik_gpd(x = data$Excess, sigma = sigma_obs_temp, xi = xi_temp))
  }
  
  # Di TIC, J adalah minus ekspektasi Hessian dari log-likelihood murni
  J <- -hessian(grad_func, param)
  
  # Hitung Score Matrix (K)
  score_matrix <- matrix(0, nrow = nrow(data), ncol = length(param))
  
  for(i in 1:nrow(data)){
    idx_stasiun <- data$Stasiun_index[i]
    
    score_i <- grad(
      function(p){
        sigma_temp <- p[1:n_stasiun]
        xi_temp <- p[length(p)]
        
        # Ambil sigma yang sesuai untuk baris observasi ini
        sigma_i <- sigma_temp[idx_stasiun]
        
        return(loglik_gpd_i(x = data$Excess[i], sigma = sigma_i, xi = xi_temp))
      },
      param
    )
    
    score_matrix[i, ] <- score_i
  }
  
  # K adalah cross-product dari score matrix
  K <- t(score_matrix) %*% score_matrix
  
  # Antisipasi error aljabar linier jika matriks J singular
  inv_J <- tryCatch(solve(J), error = function(e) {
    warning("Matriks Hessian (J) tidak bisa diinversi secara langsung. Menggunakan Pseudo-Inverse.")
    return(MASS::ginv(J))
  })
  
  TIC <- -2 * loglik + 2 * sum(diag(inv_J %*% K))
  
  return(list(
    loglik = loglik,
    TIC = TIC
  ))
}
data_exceedance <- data_exceedance %>%
  mutate(Stasiun = str_replace_all(Stasiun, "\\.", " "))
str(data_exceedance)
str(hasil_estimasi)
str(hasil_lasso)
str(hasil_ridge)
# Original
# Gabungkan data untuk mengurutkan parameter sesuai baris data_exceedance
data_tic_ori <- data_exceedance %>%
  left_join(hasil_estimasi %>% 
              dplyr::select(Stasiun, Sigma), by = "Stasiun")
# Ekstrak sigma unik (panjang = 5) sesuai urutan stasiun_unik
sigma_unik_ori <- data_tic_ori %>%
  distinct(Stasiun, .keep_all = TRUE) %>%
  pull(Sigma)
tic_ori <- hitung_TIC_Uni(
  data = data_tic_ori,
  sigma_vector = sigma_unik_ori,
  xi = unique(hasil_uji_gpd$Xi_Konstan)
)

# Fused Lasso 
data_tic_lasso <- data_exceedance %>%
  left_join(hasil_lasso %>% 
              dplyr::select(Stasiun, Sigma_FusedLasso), by = "Stasiun")
sigma_unik_lasso <- data_tic_lasso %>%
  distinct(Stasiun, .keep_all = TRUE) %>%
  pull(Sigma_FusedLasso)
tic_lasso <- hitung_TIC_Uni(
  data = data_tic_lasso,
  sigma_vector = sigma_unik_lasso,
  xi = unique(hasil_uji_gpd$Xi_Konstan)
)

# Fused Ridge 
data_tic_ridge <- data_exceedance %>%
  left_join(hasil_ridge %>% 
              dplyr::select(Stasiun, Sigma_FusedRidge), by = "Stasiun")
sigma_unik_ridge <- data_tic_ridge %>%
  distinct(Stasiun, .keep_all = TRUE) %>%
  pull(Sigma_FusedRidge)
tic_ridge <- hitung_TIC_Uni(
  data = data_tic_ridge,
  sigma_vector = sigma_unik_ridge,
  xi = unique(hasil_uji_gpd$Xi_Konstan)
)

# Perbandingan Hasil TIC untuk Model Uni
perbandingan_TIC_Uni <- data.frame(
  Metode = c("Fused Lasso", "Fused Ridge", "Original"),
  LogLik = c(tic_lasso$loglik, tic_ridge$loglik, tic_ori$loglik),
  TIC = c(tic_lasso$TIC, tic_ridge$TIC, tic_ori$TIC)
)
print(perbandingan_TIC_Uni)
write.csv(perbandingan_TIC_Uni, "D:/DontTouch!/KULEEYEEAHH SMT 7/CMEMEW/PENTING/CRHJN_TIC.csv", row.names = FALSE)





# =======================
#   VISUALISASI HEATMAP
# =======================
# Prepare Data 
data_koordinat <- data %>%
  dplyr::select(Stasiun, Latitude, Longitude) %>%
  distinct()
data_spasial <- return_level_final %>%
  dplyr::mutate(Stasiun = str_replace_all(Stasiun, "\\.", " ")) %>%
  left_join(data_koordinat, by = "Stasiun")

data_long <- data_spasial %>%
  dplyr::select(Stasiun, Longitude, Latitude, Periode_Tahun, RL_Lasso, RL_Ridge) %>%
  pivot_longer(
    cols = c(RL_Lasso, RL_Ridge),
    names_to = "Metode",
    values_to = "return_level"
  ) %>%
  mutate(
    Metode = recode(
      Metode,
      RL_Lasso = "Fused Lasso",
      RL_Ridge = "Fused Ridge"
    ),
    Periode_Tahun = factor(
      Periode_Tahun,
      levels = c(10, 25, 100),
      labels = c("10 Year", "25 Year", "100 Year")
    )
  )


# plot interpolasi
# Prepare Data 
indo <- st_read("C:/Users/ASUS/Downloads/BATAS KABUPATEN KOTA DESEMBER 2019 DUKCAPIL/BATAS KABUPATEN KOTA DESEMBER 2019 DUKCAPIL/BATAS KABUPATEN KOTA DESEMBER 2019 DUKCAPIL.shp")
indo <- st_make_valid(indo)
wilayah_pilih <- indo %>%
  filter(KAB_KOTA %in% c("TAPANULI TENGAH",
                         "PADANG LAWAS UTARA",
                         "DELI SERDANG",
                         "KOTA MEDAN",
                         "NIAS"))
grid <- st_make_grid(wilayah_pilih, cellsize = 0.05) %>%
  st_as_sf()

# Interpolasi
interpolasi_idw <- function(df){
  
  # Ambil hanya kolom yang dibutuhkan
  df_clean <- df %>%
    dplyr::select(Longitude, Latitude, return_level)
  stasiun_sf <- st_as_sf(df_clean, coords = c("Longitude", "Latitude"), crs = 4326)
  
  idw_result <- gstat::idw(
    return_level ~ 1,
    locations = stasiun_sf,
    newdata = grid
  )
  idw_sf <- st_as_sf(idw_result)
  # clipping
  idw_clip <- st_intersection(idw_sf, st_make_valid(wilayah_pilih))
  return(idw_clip)
}

idw_all <- data_long %>%
  filter(!is.na(Latitude) & !is.na(Longitude)) %>%
  group_by(Metode, Periode_Tahun) %>%
  group_split() %>%
  map_df(~{
    hasil <- interpolasi_idw(.x)
    hasil$Metode <- unique(.x$Metode)
    hasil$Periode_Tahun <- unique(.x$Periode_Tahun)
    return(hasil)
  })
idw_all <- st_make_valid(idw_all)

stasiun_sf <- st_as_sf(data_koordinat,
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
  filter(Metode == "Fused Ridge") %>%
  st_cast("POINT") %>%
  mutate(
    x = st_coordinates(.)[,1],
    y = st_coordinates(.)[,2]
  )
idw_ridge %>%
  group_by(Periode_Tahun) %>%
  summarise(min = min(var1.pred),
            max = max(var1.pred))
idwr_10 <- idw_ridge %>% filter(Periode_Tahun == "10 Year")
idwr_25 <- idw_ridge %>% filter(Periode_Tahun == "25 Year")
idwr_100 <- idw_ridge %>% filter(Periode_Tahun == "100 Year")
st_bbox(indo_sumut)

# plot df FR
plot_rl <- function(data_plot, judul){
  ggplot() +
    # Daratan Sumut
    geom_sf(data = indo_sumut,
            fill = "#cde8b6",   # hijau muda
            color = "grey70",
            linewidth = 0.2) +
    
    # Heatmap
    geom_tile(data = data_plot,
              aes(x = x, y = y, fill = var1.pred),
              width = 0.05,
              height = 0.05,
              alpha = 0.75) +
    
    # Wilayah penelitian
    geom_sf(data = wilayah_pilih,
            fill = NA,
            color = "black",
            linewidth = 0.4) +
    
    # Titik stasiun
    geom_point(data = data_spasial,
               aes(x = Longitude, y = Latitude),
               color = "red",
               size = 2) +
    geom_text(data = data_spasial,
              aes(x = Longitude, y = Latitude, label = Stasiun),
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
      title = judul,
      x = "Longitude",
      y = "Latitude"
    ) +
    
    theme_minimal() +
    theme(
      # Laut
      panel.background = element_rect(fill = "#b9e3f9", color = NA),
      # Grid
      panel.grid.major = element_line(color = "white", linewidth = 0.1),
      # Judul
      plot.title = element_text(hjust = 0.5)
    )
}


plot_10 <- plot_rl(idwr_10, "Return Level 10 Year (Fused Ridge)")
plot_10
plot_25 <- plot_rl(idwr_25, "Return Level 25 Year (Fused Ridge)")
plot_25
plot_100 <- plot_rl(idwr_100, "Return Level 100 Year (Fused Ridge)")
plot_100


idw_lasso <- idw_all %>%
  filter(Metode == "Fused Lasso") %>%
  st_cast("POINT") %>%
  mutate(
    x = st_coordinates(.)[,1],
    y = st_coordinates(.)[,2]
  )
idw_lasso %>%
  group_by(Periode_Tahun) %>%
  summarise(min = min(var1.pred),
            max = max(var1.pred))
idwl_10 <- idw_lasso %>% filter(Periode_Tahun == "10 Year")
idwl_25 <- idw_lasso %>% filter(Periode_Tahun == "25 Year")
idwl_100 <- idw_lasso %>% filter(Periode_Tahun == "100 Year")
st_bbox(indo_sumut)

plot_10 <- plot_rl(idwl_10, "Return Level 10 Year (Fused Lasso)")
plot_10
plot_25 <- plot_rl(idwl_25, "Return Level 25 Year (Fused Lasso)")
plot_25
plot_100 <- plot_rl(idwl_100, "Return Level 100 Year (Fused Lasso)")
plot_100




















