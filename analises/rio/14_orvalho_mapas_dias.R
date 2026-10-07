# ==================================================================
# 14_orvalho_mapas_dias.R
# Dew point maps for every day of the month, from a fit already saved by
# 12_orvalho_ajuste.R (nothing is refitted; takes seconds).
#
#   orvalho_dias_<tag>.png       Td of each day, all on the same colour scale
#   orvalho_anomalia_<tag>.png   each day minus the month's mean map, with the
#                                city-wide level removed: how the PATTERN changes
#   orvalho_serie_<tag>.png      city mean and spread (p05-p95) per day, with
#                                the reference stations' daily mean
# Also prints how high the stations go compared with the map (above the
# highest station the elevation effect is extrapolated).
#
# Rode de dentro de analises/rio:
#   Rscript 14_orvalho_mapas_dias.R                       (buffer 10 km)
#   Rscript 14_orvalho_mapas_dias.R tag=td_median_2026-01 (base, 40 km)
# ==================================================================

library(terra)
library(bruce)
source("config.R")

args <- commandArgs(trailingOnly = TRUE)
arg <- function(nome, padrao) {
  a <- grep(paste0("^", nome, "="), args, value = TRUE)
  if (length(a)) sub("^[^=]*=", "", a[1]) else padrao
}
TAG <- arg("tag", "td_median_2026-01_buf10km")
DIR_OUT <- file.path(RES, "orvalho")
arq <- file.path(DIR_OUT, sprintf("bruce_%s.rds", TAG))
if (!file.exists(arq)) {
  stop("Nao achei ", arq, ". Ajustes salvos: ",
       paste(list.files(DIR_OUT, "^bruce_td.*\\.rds$"), collapse = ", "), call. = FALSE)
}
res <- readRDS(arq)
td <- unwrap(res$pred_mean)
dias <- names(td)
est <- res$fit$data$stations
obs <- res$fit$data$obs
cat("Ajuste:", TAG, "-", nlyr(td), "dias\n")

# ---- 1. every day, same scale ---------------------------------------------------
lim <- quantile(values(td), c(0.005, 0.995), na.rm = TRUE)
cores <- hcl.colors(50, "YlOrRd", rev = TRUE)
nc <- 6; nl <- ceiling(nlyr(td) / nc)
png(file.path(DIR_OUT, sprintf("orvalho_dias_%s.png", TAG)), width = 2400, height = 330 * nl, res = 110)
par(mfrow = c(nl, nc), mar = c(0.5, 0.5, 1.6, 3.5))
for (i in seq_len(nlyr(td))) {
  plot(td[[i]], range = lim, col = cores, axes = FALSE, main = dias[i], legend = i %% nc == 0 | i == nlyr(td))
}
dev.off()

# ---- 2. anomaly: day minus the month's mean map -----------------------------------
# The whole city getting more or less humid (the level) is in the series plot;
# here each day's city mean is removed, so only the change of PATTERN shows.
media <- mean(td)
anom <- td - media
anom <- anom - global(anom, "mean", na.rm = TRUE)[, 1]
a <- max(abs(quantile(values(anom), c(0.005, 0.995), na.rm = TRUE)))
cores_a <- hcl.colors(50, "Blue-Red 3")
png(file.path(DIR_OUT, sprintf("orvalho_anomalia_%s.png", TAG)), width = 2400, height = 330 * (nl + 1), res = 110)
par(mfrow = c(nl + 1, nc), mar = c(0.5, 0.5, 1.6, 3.5))
plot(media, col = cores, axes = FALSE, main = "Media do mes")
for (k in 2:nc) plot.new()
for (i in seq_len(nlyr(anom))) {
  plot(anom[[i]], range = c(-a, a), col = cores_a, axes = FALSE, main = paste(dias[i], "(desenho)"),
       legend = i %% nc == 0 | i == nlyr(anom))
}
dev.off()

# ---- 3. city mean per day, and the reference stations -------------------------------
q <- t(sapply(seq_len(nlyr(td)), function(i)
  quantile(values(td[[i]]), c(0.05, 0.5, 0.95), na.rm = TRUE)))
ref <- obs[obs$source %in% c("inmet", "airport"), ]
ref_dia <- tapply(ref$value, ref$date, mean)
x <- as.Date(dias)
png(file.path(DIR_OUT, sprintf("orvalho_serie_%s.png", TAG)), width = 1300, height = 550, res = 110)
par(mar = c(4, 4.5, 2, 1))
plot(x, q[, 2], type = "n", ylim = range(c(q, ref_dia), na.rm = TRUE), xlab = "", ylab = "Td (C)",
     main = "Ponto de orvalho por dia: mapa (mediana e faixa p05-p95) e estacoes de referencia")
polygon(c(x, rev(x)), c(q[, 1], rev(q[, 3])), col = adjustcolor("orange", 0.3), border = NA)
lines(x, q[, 2], lwd = 2, col = "darkorange3")
points(as.Date(names(ref_dia)), ref_dia, pch = 19, col = "steelblue")
legend("bottomleft", c("mapa: mediana", "mapa: p05-p95", "media INMET + aeroportos"),
       col = c("darkorange3", adjustcolor("orange", 0.5), "steelblue"), lwd = c(2, 8, NA),
       pch = c(NA, NA, 19), bty = "n")
dev.off()

cat("\n=== Variacao de um dia para o outro ===\n")
cat(sprintf("Media da cidade: de %.1f a %.1f C ao longo do mes\n", min(q[, 2]), max(q[, 2])))
dp_anom <- global(anom, "sd", na.rm = TRUE)[, 1]
cat(sprintf("Desvio espacial da anomalia (o quanto o DESENHO muda no dia): mediana %.2f C, maximo %.2f C (%s)\n",
            median(dp_anom), max(dp_anom), dias[which.max(dp_anom)]))

# ---- 4. how high the stations go ----------------------------------------------------
arq_dem <- file.path(DIR_OUT, "cache", "dem.tif")
if (file.exists(arq_dem)) {
  dem <- rast(arq_dem)
  z_est <- extract(dem, vect(est, geom = c("x", "y"), crs = crs(dem)))[, 2]
  z_map <- extract(dem, as.points(td[[1]]))[, 2]
  cat("\n=== Altitude: estacoes x mapa ===\n")
  cat(sprintf("Estacao mais alta: %.0f m (%s)\n", max(z_est, na.rm = TRUE), est$station_id[which.max(z_est)]))
  cat(sprintf("Celulas do mapa acima da estacao mais alta: %.1f%% (ate %.0f m)\n",
              100 * mean(z_map > max(z_est, na.rm = TRUE), na.rm = TRUE), max(z_map, na.rm = TRUE)))
  cat("Nessas celulas o efeito da elevacao e extrapolado (reta alem dos dados).\n")
}
cat("\nFiguras em", DIR_OUT, ": orvalho_dias_, orvalho_anomalia_, orvalho_serie_", TAG, "\n")
