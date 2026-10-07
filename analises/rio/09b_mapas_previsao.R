# ==================================================================
# 09b_mapas_previsao.R
# Maps of the next-day forecasts of 09_previsao_dia_seguinte.R.
#
# For each window: the last day of the window (with data), the forecast
# for the next day (no data) and its uncertainty, on the article's 500 m
# prediction grid (dp in each package, covariates already standardised).
# The forecast panel shows the stations with what they actually measured
# that day, on the same colour scale (station readings include each
# station's own offset; the map is the reference-scale surface).
#
# Refits each window as in 09 (theta fixed, ~1 min each).
# Rode de dentro de analises/rio, depois do 09:
#   Rscript 09b_mapas_previsao.R > log_mapas_previsao.log 2>&1
# ==================================================================

library(bruce)
library(sf)
library(terra)
suppressMessages(library(INLA))
source("config.R")
inla.setOption(fmesher.evolution.warn = FALSE)

VAR <- "median"
PACOTES <- path.expand("~/MEGA/CIE things/pacotes_artigo")
DIR_THETA <- file.path(RES, "theta_emprestado")
DIR_PREV <- file.path(RES, "previsao_dia_seguinte")
DIR <- file.path(DIR_PREV, "mapas")
dir.create(DIR, showWarnings = FALSE, recursive = TRUE)
N_AMOSTRAS <- 100

arqs <- list.files(PACOTES, sprintf("^ARTIGO_pacote_%s_[0-9]{4}-[0-9]{2}\\.rds$", VAR), full.names = TRUE)
janelas <- sub(".*_([0-9]{4}-[0-9]{2})\\.rds$", "\\1", basename(arqs))
names(arqs) <- janelas
manifesto <- read.csv(file.path(PACOTES, "MANIFESTO_pacotes_artigo.csv"))
manifesto <- manifesto[manifesto$stat == VAR, ]
spec <- bruce_spec(range_prior = 15000, intercept_prior = c(0, 0.001), on_collapse = "warn",
                   int_strategy = "eb")
spec_janela <- function(j) {
  s <- spec
  m <- manifesto[manifesto$mes == j, ]
  if (nrow(m) == 1 && m$grouprho == "fixo") s$group_rho_fixed <- m$grouprho_valor
  s
}

# The article's inputs for one window (same as 04, 08, 09, 10)
entradas <- function(orig) {
  est <- st_drop_geometry(orig$est_all)
  src <- unname(c(ofic = "official", aura = "aurassure", wunder = "wunder")[est$origem])
  src[est$estacao_id %in% orig$ESTACOES_AEROPORTO] <- "airport"
  stations <- data.frame(station_id = est$estacao_id, source = src, x = est$x, y = est$y)
  stations <- stations[!duplicated(stations$station_id), ]
  obs <- data.frame(station_id = est$estacao_id, date = as.Date(est$data),
                    value = est[[orig$col_estacao]])
  bg <- st_drop_geometry(orig$d_all)
  background <- data.frame(x = bg$x, y = bg$y, date = as.Date(bg$date), value = bg[[orig$var_era5]])
  sources <- data.frame(source = c("official", "airport", "aurassure", "wunder"),
                        bias = c(FALSE, TRUE, TRUE, TRUE), obs_sd = c(0.3, 0.3, NA, NA))
  sources <- sources[sources$source %in% stations$source, ]
  d <- bruce_data(obs, stations, sources, crs = 31983, background = background, background_sd = 0.7)
  st_cov <- data.frame(station_id = est$estacao_id, elev = est$elev_std, urb = est$urb_std,
                       svf = est$svf_std)
  st_cov <- st_cov[!duplicated(st_cov$station_id), ]
  st_cov$svf[is.na(st_cov$svf)] <- 0
  bg_cov <- data.frame(x = bg$x, y = bg$y, elev = bg$elev_std, urb = bg$urb_std, svf = 0)
  bg_cov <- bg_cov[!duplicated(bg_cov[, c("x", "y")]), ]
  d$covariates <- list(names = c("elev", "urb", "svf"), rasters = NULL, scaling = NULL,
                       station = st_cov, background = bg_cov)
  d
}

# The article's prediction grid: x, y and standardised covariates
grade_artigo <- function(orig) {
  dp <- as.data.frame(orig$dp)
  if (inherits(orig$dp, "sf")) dp <- sf::st_drop_geometry(orig$dp)
  precisa <- c("x", "y", "elev_std", "urb_std", "svf_std")
  falta <- setdiff(precisa, names(dp))
  if (length(falta)) stop("A grade dp nao tem: ", paste(falta, collapse = ", "),
                          ". Colunas disponiveis: ", paste(names(dp), collapse = ", "))
  dp <- unique(dp[, precisa])
  dp$svf_std[is.na(dp$svf_std)] <- 0
  dp
}

# Mean and sd of the reference-scale surface on the grid for day g
superficie_grade <- function(fit, dp, g) {
  X <- cbind(intercept = 1, elev = dp$elev_std, urb = dp$urb_std, svf = dp$svf_std)[, c("intercept", fit$cov_names)]
  A <- INLA::inla.spde.make.A(fit$mesh, loc = cbind(dp$x, dp$y))
  comps <- bruce:::latent_components(fit, N_AMOSTRAS, NULL)
  P <- bruce:::eval_predictor(fit, comps, X, A, g)
  data.frame(x = dp$x, y = dp$y, media = rowMeans(P), sd = apply(P, 1, stats::sd))
}
para_raster <- function(df, col) terra::rast(df[, c("x", "y", col)], type = "xyz", crs = "EPSG:31983")

paineis <- list()
for (j in janelas) {
  f_prev <- file.path(DIR_PREV, sprintf("previsao_%s.rds", j))
  f_theta <- file.path(DIR_THETA, sprintf("proprio_%s.rds", j))
  if (!file.exists(f_prev) || !file.exists(f_theta)) { cat(sprintf("[%s] falta o 09 ou o 08, pulado\n", j)); next }
  prev <- readRDS(f_prev)
  orig <- readRDS(arqs[[j]])
  d <- entradas(orig)
  dp <- grade_artigo(orig)
  ultimo <- max(d$obs$date); alvo <- ultimo + 1
  cat(sprintf("[%s] mapas de %s (ultimo dia) e %s (previsao)\n", j, format(ultimo), format(alvo)))

  fit <- bruce_fit(d, orig$mesh, spec_janela(j), window = c(min(d$obs$date), alvo),
                   theta_init = readRDS(f_theta)$theta, fix_theta = TRUE, keep_config = TRUE)
  s_ult <- superficie_grade(fit, dp, which(fit$days == ultimo))
  s_prv <- superficie_grade(fit, dp, which(fit$days == alvo))
  r <- c(para_raster(s_ult, "media"), para_raster(s_prv, "media"), para_raster(s_prv, "sd"))
  names(r) <- c(paste0("ultimo_", ultimo), paste0("previsao_", alvo), paste0("sd_previsao_", alvo))
  writeRaster(r, file.path(DIR, sprintf("mapas_%s.tif", j)), overwrite = TRUE)

  xy <- unique(d$obs[, c("station_id", "x", "y")])
  prev$x <- xy$x[match(prev$station_id, xy$station_id)]
  prev$y <- xy$y[match(prev$station_id, xy$station_id)]
  paineis[[j]] <- list(r = r, prev = prev, ultimo = ultimo, alvo = alvo)

  lim <- range(c(values(r[[1:2]]), prev$obs), na.rm = TRUE)
  pal <- hcl.colors(60, "YlOrRd", rev = TRUE)
  cor_de <- function(v) pal[pmax(1, pmin(60, round(1 + 59 * (v - lim[1]) / diff(lim))))]
  png(file.path(DIR, sprintf("mapas_%s.png", j)), width = 1650, height = 560, res = 110)
  par(mfrow = c(1, 3), mar = c(1, 1, 3, 4))
  plot(r[[1]], col = pal, range = lim, axes = FALSE, main = sprintf("Ultimo dia com dados (%s)", format(ultimo)))
  plot(r[[2]], col = pal, range = lim, axes = FALSE,
       main = sprintf("Previsao para %s\n(pontos: medido nas estacoes)", format(alvo)))
  points(prev$x, prev$y, pch = 21, bg = cor_de(prev$obs), cex = 1.4)
  plot(r[[3]], col = hcl.colors(60, "Purples", rev = TRUE), axes = FALSE,
       main = sprintf("Incerteza da previsao (sd, C)"))
  dev.off()
  rm(fit, orig, d); gc()
}

# All forecasts side by side, each with its own colour scale
if (length(paineis)) {
  nc <- min(3, length(paineis)); nl <- ceiling(length(paineis) / nc)
  png(file.path(DIR, "previsoes_todas.png"), width = 600 * nc, height = 520 * nl, res = 110)
  par(mfrow = c(nl, nc), mar = c(1, 1, 3, 4))
  for (j in names(paineis)) {
    p <- paineis[[j]]
    lim <- range(c(values(p$r[[2]]), p$prev$obs), na.rm = TRUE)
    pal <- hcl.colors(60, "YlOrRd", rev = TRUE)
    plot(p$r[[2]], col = pal, range = lim, axes = FALSE,
         main = sprintf("%s (RMSE %.2f C)", format(p$alvo), sqrt(mean((p$prev$prev - p$prev$obs)^2))))
    points(p$prev$x, p$prev$y, pch = 21, cex = 1.2,
           bg = pal[pmax(1, pmin(60, round(1 + 59 * (p$prev$obs - lim[1]) / diff(lim))))])
  }
  dev.off()
}
cat("Mapas em", DIR, "\n")
