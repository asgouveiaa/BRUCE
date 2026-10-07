# ==================================================================
# 16_indice_calor.R
# Heat index map for one day: daily MAXIMUM temperature (article's model)
# + daily median DEW POINT (12_orvalho_ajuste.R).
#
# Why max temperature with the day's dew point: the dew point barely
# changes during the day, so Tmax with the day's dew point gives the
# humidity of the hottest hour, the one that matters for heat stress.
#
#   1. Tmax: refits the article's Tmax model for the target day +- 7 days
#      (inputs from ARTIGO_pacote_max_<month>.rds, eb, ~ a few minutes;
#      saved, so a rerun skips it).
#   2. Dew point: the saved fit from 12 (nothing refitted).
#   3. 200 joint draws of each map; for each pair of draws:
#      relative humidity (Magnus) and the heat index (US NWS formula,
#      Rothfusz with its adjustments). Mean, sd and 5-95% of the heat index.
#   4. Check at the official stations: heat index from what they measured
#      (Tmax and dew point) against the map.
#
# The two models are fitted separately, so their draws are independent;
# real errors of T and Td at one place may be correlated (the interval is
# approximate, the mean map is not affected). Where a draw gives Td > T,
# relative humidity is capped at 100%.
#
# Rode de dentro de analises/rio:
#   Rscript 16_indice_calor.R > log_indice_calor.log 2>&1
#   Rscript 16_indice_calor.R dia=2023-11-18 td=td_median_2023-11_buf10km_grho0.9_sem2
# ==================================================================

library(bruce)
library(sf)
library(terra)
suppressMessages(library(INLA))
source("config.R")
inla.setOption(fmesher.evolution.warn = FALSE)

args <- commandArgs(trailingOnly = TRUE)
arg <- function(nome, padrao) {
  a <- grep(paste0("^", nome, "="), args, value = TRUE)
  if (length(a)) sub("^[^=]*=", "", a[1]) else padrao
}
DIA <- as.Date(arg("dia", "2023-11-18"))
TAG_TD <- arg("td", "td_median_2023-11_buf10km_grho0.9_sem2")
MEIA_JANELA <- 7                 # Tmax fitted on DIA +- 7 days
N <- 200                         # draws
MES <- format(DIA, "%Y-%m")
PACOTES <- path.expand("~/MEGA/CIE things/pacotes_artigo")
DIR_TD <- file.path(RES, "orvalho")
DIR <- file.path(RES, "indice_calor")
dir.create(DIR, showWarnings = FALSE, recursive = TRUE)
cat("Dia:", format(DIA), " ponto de orvalho:", TAG_TD, "\n")

# ---- heat index (US NWS) ------------------------------------------------------------
umidade_relativa <- function(t, td) {
  e <- function(x) 6.112 * exp(17.62 * x / (243.12 + x))
  pmin(100 * e(td) / e(t), 100)        # first argument first: keeps the matrix shape
}
indice_calor <- function(t, ur) {            # t in C, ur in %; returns C
  tf <- t * 9 / 5 + 32
  simples <- 0.5 * (tf + 61 + (tf - 68) * 1.2 + ur * 0.094)
  hi <- -42.379 + 2.04901523 * tf + 10.14333127 * ur - 0.22475541 * tf * ur -
    0.00683783 * tf^2 - 0.05481717 * ur^2 + 0.00122874 * tf^2 * ur +
    0.00085282 * tf * ur^2 - 0.00000199 * tf^2 * ur^2
  aj1 <- ur < 13 & tf >= 80 & tf <= 112
  hi[aj1] <- hi[aj1] - ((13 - ur[aj1]) / 4) * sqrt((17 - abs(tf[aj1] - 95)) / 17)
  aj2 <- ur > 85 & tf >= 80 & tf <= 87
  hi[aj2] <- hi[aj2] + ((ur[aj2] - 85) / 10) * ((87 - tf[aj2]) / 5)
  usa_simples <- (simples + tf) / 2 < 80
  hi[usa_simples] <- simples[usa_simples]
  (hi - 32) * 5 / 9
}

# ---- 1. Tmax: the article's model ------------------------------------------------------
arq_pac <- file.path(PACOTES, sprintf("ARTIGO_pacote_max_%s.rds", MES))
if (!file.exists(arq_pac)) stop("Nao achei ", arq_pac, call. = FALSE)
orig <- readRDS(arq_pac)

entradas <- function(orig) {                # same as 04, 08, 09, 10
  est <- st_drop_geometry(orig$est_all)
  src <- unname(c(ofic = "official", aura = "aurassure", wunder = "wunder")[est$origem])
  src[est$estacao_id %in% orig$ESTACOES_AEROPORTO] <- "airport"
  stations <- data.frame(station_id = est$estacao_id, source = src, x = est$x, y = est$y)
  stations <- stations[!duplicated(stations$station_id), ]
  obs <- data.frame(station_id = est$estacao_id, date = as.Date(est$data), value = est[[orig$col_estacao]])
  bg <- st_drop_geometry(orig$d_all)
  background <- data.frame(x = bg$x, y = bg$y, date = as.Date(bg$date), value = bg[[orig$var_era5]])
  sources <- data.frame(source = c("official", "airport", "aurassure", "wunder"),
                        bias = c(FALSE, TRUE, TRUE, TRUE), obs_sd = c(0.3, 0.3, NA, NA))
  sources <- sources[sources$source %in% stations$source, ]
  d <- bruce_data(obs, stations, sources, crs = 31983, background = background, background_sd = 0.7)
  st_cov <- data.frame(station_id = est$estacao_id, elev = est$elev_std, urb = est$urb_std, svf = est$svf_std)
  st_cov <- st_cov[!duplicated(st_cov$station_id), ]
  st_cov$svf[is.na(st_cov$svf)] <- 0
  bg_cov <- data.frame(x = bg$x, y = bg$y, elev = bg$elev_std, urb = bg$urb_std, svf = 0)
  bg_cov <- bg_cov[!duplicated(bg_cov[, c("x", "y")]), ]
  d$covariates <- list(names = c("elev", "urb", "svf"), rasters = NULL, scaling = NULL,
                       station = st_cov, background = bg_cov)
  d
}
d_t <- entradas(orig)
dp <- as.data.frame(if (inherits(orig$dp, "sf")) st_drop_geometry(orig$dp) else orig$dp)
dp <- unique(dp[, c("x", "y", "elev_std", "urb_std", "svf_std")])
dp$svf_std[is.na(dp$svf_std)] <- 0

arq_fit_t <- file.path(DIR, sprintf("fit_tmax_%s.rds", format(DIA)))
if (file.exists(arq_fit_t)) {
  fit_t <- readRDS(arq_fit_t); cat("Tmax: ajuste salvo carregado\n")
} else {
  spec_t <- bruce_spec(range_prior = 15000, intercept_prior = c(0, 0.001), int_strategy = "eb")
  man <- file.path(PACOTES, "MANIFESTO_pacotes_artigo.csv")
  if (file.exists(man)) {
    m <- read.csv(man); m <- m[m$stat == "max" & m$mes == MES, ]
    if (nrow(m) == 1 && m$grouprho == "fixo") spec_t$group_rho_fixed <- m$grouprho_valor
  }
  cat("Tmax: ajustando", format(DIA - MEIA_JANELA), "a", format(DIA + MEIA_JANELA), "...\n")
  fit_t <- bruce_fit(d_t, orig$mesh, spec_t, window = c(DIA - MEIA_JANELA, DIA + MEIA_JANELA),
                     keep_config = TRUE)
  saveRDS(fit_t, arq_fit_t)
}
invisible(bruce_check(fit_t))
if (!DIA %in% fit_t$days) stop("O dia ", DIA, " nao esta no ajuste da Tmax", call. = FALSE)
X_t <- cbind(intercept = 1, elev = dp$elev_std, urb = dp$urb_std, svf = dp$svf_std)[, c("intercept", fit_t$cov_names)]
A_t <- INLA::inla.spde.make.A(fit_t$mesh, loc = cbind(dp$x, dp$y))
P_t <- bruce:::eval_predictor(fit_t, bruce:::latent_components(fit_t, N, NULL), X_t, A_t,
                              which(fit_t$days == DIA))

# ---- 2. dew point on the same grid ------------------------------------------------------
arq_td <- file.path(DIR_TD, sprintf("bruce_%s.rds", TAG_TD))
if (!file.exists(arq_td)) stop("Nao achei ", arq_td, call. = FALSE)
fit_td <- readRDS(arq_td)$fit
grade <- rast(cbind(dp[, c("x", "y")], v = 1), type = "xyz", crs = "EPSG:31983")
pred_td <- bruce_predict(fit_td, grid = grade, dates = DIA, n_samples = N, keep_samples = TRUE,
                         resample = c(elev = "average", urb = "near", dist_mar = "average"), seed = NULL)
P_td_cel <- pred_td$samples[[1]]                       # rows = pred_td$cells
cel_dp <- cellFromXY(grade, as.matrix(dp[, c("x", "y")]))
P_td <- P_td_cel[match(cel_dp, pred_td$cells), , drop = FALSE]
ok <- stats::complete.cases(P_td) & stats::complete.cases(P_t)
cat(sprintf("Celulas: %d (com Tmax e ponto de orvalho: %d)\n", nrow(dp), sum(ok)))

# ---- 3. heat index from each pair of draws ---------------------------------------------
UR <- umidade_relativa(P_t, P_td)
HI <- indice_calor(P_t, UR)
res <- data.frame(x = dp$x, y = dp$y,
                  tmax = rowMeans(P_t), td = rowMeans(P_td), ur = rowMeans(UR),
                  ic = rowMeans(HI), ic_sd = apply(HI, 1, sd),
                  ic_p05 = apply(HI, 1, quantile, 0.05, na.rm = TRUE),
                  ic_p95 = apply(HI, 1, quantile, 0.95, na.rm = TRUE))
res <- res[ok, ]
r <- rast(res, type = "xyz", crs = "EPSG:31983")
writeRaster(r, file.path(DIR, sprintf("indice_calor_%s.tif", format(DIA))), overwrite = TRUE)

cat(sprintf("\n=== %s: resumo sobre a cidade (media das celulas; min - max) ===\n", format(DIA)))
for (v in c("tmax", "td", "ur", "ic", "ic_sd")) {
  cat(sprintf("  %-6s %6.1f   (%5.1f - %5.1f)\n", v, mean(res[[v]]), min(res[[v]]), max(res[[v]])))
}
faixas <- cut(res$ic, c(-Inf, 27, 32, 41, 54, Inf),
              c("< 27 normal", "27-32 cautela", "32-41 cautela extrema", "41-54 perigo", ">= 54 perigo extremo"))
cat("\nArea da cidade por faixa do indice de calor (NWS):\n")
print(round(100 * prop.table(table(faixas)), 1))

# ---- 4. check at the official stations ----------------------------------------------
tryCatch({
  tobs <- d_t$obs[d_t$obs$date == DIA & d_t$obs$source %in% c("official", "airport"), ]
  tdc <- readRDS(file.path(DIR_TD, "cache", sprintf("td_obs_median_%s.rds", MES)))
  tdc <- tdc[as.Date(tdc$data) == DIA & tdc$origem %in% c("inmet", "redemet") & !is.na(tdc$median_orvalho), ]
  xy_td <- st_coordinates(st_transform(st_as_sf(tdc, coords = c("longitude", "latitude"), crs = 4326), 31983))
  linhas <- lapply(seq_len(nrow(tobs)), function(i) {
    dd <- sqrt((xy_td[, 1] - tobs$x[i])^2 + (xy_td[, 2] - tobs$y[i])^2)
    if (!length(dd) || min(dd) > 2000) return(NULL)
    j <- which.min(dd)
    k <- which.min((res$x - tobs$x[i])^2 + (res$y - tobs$y[i])^2)
    ic_obs <- indice_calor(tobs$value[i], umidade_relativa(tobs$value[i], tdc$median_orvalho[j]))
    data.frame(estacao = tobs$station_id[i], tmax_obs = tobs$value[i], td_obs = tdc$median_orvalho[j],
               ic_obs = ic_obs, ic_mapa = res$ic[k], ic_p05 = res$ic_p05[k], ic_p95 = res$ic_p95[k])
  })
  tab <- do.call(rbind, linhas)
  if (!is.null(tab)) {
    num <- vapply(tab, is.numeric, TRUE); tab[num] <- lapply(tab[num], round, 1)
    cat("\n=== Estacoes oficiais: indice de calor do que mediram x mapa ===\n")
    print(tab, row.names = FALSE)
    cat("(o mapa e a superficie de referencia; a leitura de cada estacao tem seu proprio desvio)\n")
  }
}, error = function(e) cat("\n(checagem nas estacoes nao rodou: ", conditionMessage(e), ")\n"))

# ---- 5. figure ----------------------------------------------------------------------
png(file.path(DIR, sprintf("indice_calor_%s.png", format(DIA))), width = 2000, height = 1050, res = 110)
par(mfrow = c(2, 3), mar = c(1, 1, 3, 5))
plot(r$tmax, col = hcl.colors(50, "YlOrRd", rev = TRUE), axes = FALSE, main = "Temperatura maxima (C)")
plot(r$td, col = hcl.colors(50, "YlGnBu", rev = TRUE), axes = FALSE, main = "Ponto de orvalho (C)")
plot(r$ur, col = hcl.colors(50, "Blues", rev = TRUE), axes = FALSE, main = "Umidade relativa na hora mais quente (%)")
plot(r$ic, col = hcl.colors(50, "Inferno", rev = TRUE), axes = FALSE,
     main = sprintf("Indice de calor (C) - %s", format(DIA)))
# On an extreme day the whole city falls in one NWS class, so this panel shows
# where it is hotter or cooler than the city as a whole (the NWS classes are
# in the log)
dif <- r$ic - mean(res$ic)
lim_d <- max(abs(range(values(dif), na.rm = TRUE)))
plot(dif, col = hcl.colors(50, "Blue-Red 3"), range = c(-lim_d, lim_d), axes = FALSE,
     main = sprintf("Indice de calor menos a media da cidade (%.1f C)", mean(res$ic)))
plot(r$ic_sd, col = hcl.colors(50, "Purples", rev = TRUE), axes = FALSE, main = "Incerteza do indice (sd, C)")
dev.off()
cat("\nSalvos em", DIR, ":", sprintf("indice_calor_%s.tif e .png", format(DIA)), "\n")
