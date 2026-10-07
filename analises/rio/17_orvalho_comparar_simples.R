# ==================================================================
# 17_orvalho_comparar_simples.R
# Is bruce worth it for dew point? The same leave-station-out test as bruce
# (same official stations held out, same days), with simple methods:
#
#   era5            ERA5-Land d2m of the nearest cell, as it is
#   era5_vies       + the day's mean difference at the other official
#                   stations (removes ERA5's dryness, day by day)
#   era5_alt        ERA5 + regression on all the other stations:
#                   day + network + elevation difference (station - ERA5 cell)
#   idw             inverse-distance weighting of the other stations'
#                   readings, citizen networks corrected by their offset
#   era5_alt_idw    era5_alt + IDW of its station residuals (the simple
#                   version of "covariates + spatial field")
#   bruce           the LSO already saved by 12_orvalho_ajuste.R
#
# Everything the simple methods use comes from the fit saved by 12 (same
# stations, same ERA5 cells, same cleaning). Fits nothing; takes seconds.
#
# If bruce wins clearly, the spatial model earns its place for dew point.
# If a simple method ties, dew point can be mapped the simple way and bruce
# stays for temperature. Both are useful answers.
#
# Rode de dentro de analises/rio:
#   Rscript 17_orvalho_comparar_simples.R
#   Rscript 17_orvalho_comparar_simples.R tag=td_median_2026-01_buf10km
# ==================================================================

library(terra)
source("config.R")

args <- commandArgs(trailingOnly = TRUE)
arg <- function(nome, padrao) {
  a <- grep(paste0("^", nome, "="), args, value = TRUE)
  if (length(a)) sub("^[^=]*=", "", a[1]) else padrao
}
TAG <- arg("tag", "td_median_2023-11_buf10km_grho0.9_sem2")
DIR_OUT <- file.path(RES, "orvalho")
REF <- c("inmet", "airport")         # reference networks
POT_IDW <- 2

fit <- readRDS(file.path(DIR_OUT, sprintf("bruce_%s.rds", TAG)))$fit
arq_cv <- file.path(DIR_OUT, sprintf("bruce_lso_%s.rds", TAG))
if (!file.exists(arq_cv)) stop("Nao achei o LSO do bruce: ", arq_cv, call. = FALSE)
cv <- readRDS(arq_cv)
pb <- as.data.frame(cv$predictions)[, c("station_id", "date", "obs", "pred")]
pb$date <- as.Date(pb$date)

obs <- as.data.frame(fit$data$obs)[, c("station_id", "date", "value", "source", "x", "y")]
obs <- obs[obs$date %in% fit$days, ]
bg <- as.data.frame(fit$data$background)
bg <- bg[bg$date %in% fit$days, ]
cel <- unique(bg[, c("x", "y")])
cat(sprintf("Ajuste %s: %d estacoes, %d dias, %d celulas ERA5\n", TAG,
            length(unique(obs$station_id)), length(fit$days), nrow(cel)))

# elevation of each station and of each ERA5 cell (mean within 5 km, about
# the size of an ERA5-Land cell)
dem <- rast(file.path(DIR_OUT, "cache", "dem.tif"))
est <- unique(obs[, c("station_id", "source", "x", "y")])
est$elev <- extract(dem, vect(est, geom = c("x", "y"), crs = crs(dem)))[, 2]
cel$elev <- vapply(seq_len(nrow(cel)), function(i) {
  b <- buffer(vect(cel[i, ], geom = c("x", "y"), crs = crs(dem)), 5000)
  mean(extract(dem, b)[, 2], na.rm = TRUE)
}, 0)
est$cel <- vapply(seq_len(nrow(est)), function(i)
  which.min((cel$x - est$x[i])^2 + (cel$y - est$y[i])^2), 1L)
est$delta_elev <- est$elev - cel$elev[est$cel]

bg$cel <- match(paste(bg$x, bg$y), paste(cel$x, cel$y))
obs <- merge(obs, est[, c("station_id", "cel", "delta_elev")], by = "station_id")
obs$era5 <- bg$value[match(paste(obs$cel, obs$date), paste(bg$cel, bg$date))]
obs <- obs[!is.na(obs$era5), ]
obs$dif <- obs$value - obs$era5

idw <- function(x0, y0, x, y, v) {
  w <- 1 / pmax(sqrt((x - x0)^2 + (y - y0)^2), 100)^POT_IDW
  sum(w * v) / sum(w)
}

# ---- leave each held-out station of bruce's LSO out --------------------------------------
fora <- unique(pb$station_id)
linhas <- list()
for (s in fora) {
  tr <- obs[obs$station_id != s, ]
  te <- obs[obs$station_id == s, ]
  te <- te[te$date %in% pb$date[pb$station_id == s], ]
  if (!nrow(te)) next
  # day + network + elevation difference, on all the other stations (reference = 0)
  tr$rede <- ifelse(tr$source %in% REF, "ref", tr$source)
  tr$rede <- relevel(factor(tr$rede), "ref")
  m <- lm(dif ~ factor(date) + rede + delta_elev, data = tr)
  off_rede <- c(ref = 0, coef(m)[grep("^rede", names(coef(m)))])
  names(off_rede) <- sub("^rede", "", names(off_rede))
  tr$corr <- tr$value - off_rede[ifelse(tr$source %in% REF, "ref", tr$source)]
  tr$res_alt <- residuals(m)
  for (k in seq_len(nrow(te))) {
    dd <- te$date[k]
    ref_d <- tr[tr$date == dd & tr$source %in% REF, ]
    all_d <- tr[tr$date == dd, ]
    p_alt <- tryCatch(predict(m, newdata = data.frame(date = dd, rede = factor("ref", levels = levels(tr$rede)),
                                                      delta_elev = te$delta_elev[k])), error = function(e) NA)
    linhas[[length(linhas) + 1]] <- data.frame(
      station_id = s, date = dd, obs = te$value[k],
      era5 = te$era5[k],
      era5_vies = te$era5[k] + if (nrow(ref_d)) mean(ref_d$dif) else NA,
      era5_alt = te$era5[k] + p_alt,
      idw = if (nrow(all_d)) idw(te$x[k], te$y[k], all_d$x, all_d$y, all_d$corr) else NA,
      era5_alt_idw = te$era5[k] + p_alt +
        if (nrow(all_d)) idw(te$x[k], te$y[k], all_d$x, all_d$y, all_d$res_alt) else 0)
  }
}
pr <- do.call(rbind, linhas)
pr <- merge(pr, pb[, c("station_id", "date", "pred")], by = c("station_id", "date"))
names(pr)[names(pr) == "pred"] <- "bruce"
metodos <- c("era5", "era5_vies", "era5_alt", "idw", "era5_alt_idw", "bruce")
pr <- pr[stats::complete.cases(pr[, c("obs", metodos)]), ]

nota <- function(p) {
  do.call(rbind, lapply(metodos, function(mt) {
    e <- p[[mt]] - p$obs
    data.frame(metodo = mt, n = length(e), rmse = round(sqrt(mean(e^2)), 3),
               mae = round(mean(abs(e)), 3), vies = round(mean(e), 3))
  }))
}
cat(sprintf("\n=== %s: estacoes oficiais deixadas de fora, mesmos dias para todos ===\n", TAG))
geral <- nota(pr)
print(geral[order(geral$rmse), ], row.names = FALSE)

cat("\nRMSE por estacao:\n")
por <- do.call(rbind, lapply(split(pr, pr$station_id), function(p) {
  v <- sapply(metodos, function(mt) round(sqrt(mean((p[[mt]] - p$obs)^2)), 2))
  data.frame(estacao = p$station_id[1], n = nrow(p), t(v))
}))
print(por, row.names = FALSE)

# does bruce beat the best simple method on more days than chance?
melhor <- geral$metodo[geral$metodo != "bruce"][which.min(geral$rmse[geral$metodo != "bruce"])]
d_err <- abs(pr$bruce - pr$obs) - abs(pr[[melhor]] - pr$obs)
tt <- t.test(d_err)
cat(sprintf("\nbruce x melhor metodo simples (%s): erro absoluto do bruce menor em %.0f%% dos casos;\n",
            melhor, 100 * mean(d_err < 0)))
cat(sprintf("diferenca media do erro absoluto %.3f C (IC95%% %.3f a %.3f; negativo = bruce melhor)\n",
            mean(d_err), tt$conf.int[1], tt$conf.int[2]))
cat("(os dias de uma mesma estacao nao sao independentes: o IC e otimista)\n")

write.csv(pr, file.path(DIR_OUT, sprintf("comparar_simples_%s.csv", TAG)), row.names = FALSE)
cat("\nPrevisoes salvas em", file.path(DIR_OUT, sprintf("comparar_simples_%s.csv", TAG)), "\n")
