# ==================================================================
# 19_rodar_mes.R
# bruce for any month of Rio de Janeiro, from start to finish: data, fit,
# checks, maps, validation and figures. Same model as the article (01) by
# default; everything goes to resultados/mes_<var>_<yyyy-mm>/.
#
#   1. Stations (query_meteorj_periodo.R) and ERA5-Land t2m, covariates
#      (elevation, urban fraction, sky view factor; cached after the first run).
#   2. Diagnostic: ERA5 minus each station, by network; the stations farthest
#      from their own network are listed (qc=3 drops citizen stations more
#      than 3 C away from their network, or with fewer than 5 days).
#   3. Fit (eb by default; estrategia=ccd as in the article) and bruce_check().
#   4. Maps of every day (mean and uncertainty), decomposition of the
#      hottest day, city series against the official stations.
#   5. Leave-station-out validation on official stations (starts from the
#      fit's theta), with RMSE, bias, CRPS and coverage, per station and in
#      a figure. lro=1 also runs leave-region-out (hours).
#
# Rode de dentro de analises/rio:
#   Rscript 19_rodar_mes.R ano=2024 mes=2 > log_mes_2024-02.log 2>&1 &
#   Rscript 19_rodar_mes.R ano=2024 mes=2 var=max qc=3 dobras=0
# Options:
#   ano=, mes=        month (default 2025-01)
#   var=median        "median", "max" or "min"
#   estrategia=eb     "eb" (faster, same results in 08) or "ccd" (article)
#   qc=3              quality control of citizen stations (default: off)
#   dobras=10         official stations held out in the LSO (0 = all; each
#                     takes ~2-6 min); cv=0 skips the validation
#   lro=1             also leave-region-out by planning area (slow)
#   resolucao=500     prediction grid (m)
# Long runs: use nohup (or tmux), so the run survives closing the terminal:
#   nohup Rscript 19_rodar_mes.R ano=2024 mes=2 > log_mes_2024-02.log 2>&1 &
# Time: data ~5 min (first run longer), fit 10-30 min, LSO ~2-6 min per
# station. Each LSO fold is saved; rerunning resumes.
# ==================================================================

library(bruce)
library(terra)
library(sf)
source("config.R")
source("query_meteorj_periodo.R")
source(ERA5_FUNCOES)
if (packageVersion("bruce") < "0.1.0.9") stop("Precisa do bruce >= 0.1.0.9")

args <- commandArgs(trailingOnly = TRUE)
arg <- function(nome, padrao) {
  a <- grep(paste0("^", nome, "="), args, value = TRUE)
  if (length(a)) sub("^[^=]*=", "", a[1]) else padrao
}
ANO <- as.integer(arg("ano", 2025)); MES <- as.integer(arg("mes", 1))
VAR <- arg("var", "median")
ESTRATEGIA <- arg("estrategia", "eb")
QC <- as.numeric(arg("qc", NA))
DOBRAS <- as.integer(arg("dobras", 10))
RODAR_CV <- arg("cv", "1") == "1"
RODAR_LRO <- arg("lro", "0") == "1"
RESOLUCAO <- as.numeric(arg("resolucao", 500))
N_AMOSTRAS <- 100
CRS_UTM <- 31983
AIRPORTS <- c("SBJR", "SBGL", "SBAF", "SBSC")

tag <- sprintf("%s_%d-%02d", VAR, ANO, MES)
DIR <- file.path(RES, sprintf("mes_%s", tag))
dir.create(DIR, showWarnings = FALSE, recursive = TRUE)
DIR_COV <- file.path(RES, "cache_covariaveis")
dir.create(DIR_COV, showWarnings = FALSE)
col_est <- paste0(VAR, "_temperatura")
var_era5 <- paste0("t2m_", VAR)
data_ini <- as.Date(sprintf("%d-%02d-01", ANO, MES))
data_fim <- seq(data_ini, by = "month", length.out = 2)[2] - 1
cat(sprintf("== bruce %s | %s a %s | %s%s ==\n", VAR, data_ini, data_fim, ESTRATEGIA,
            if (!is.na(QC)) sprintf(" | qc=%g", QC) else ""))

mrj <- geobr::read_municipality(code_muni = 3304557, year = 2020, showProgress = FALSE) |> st_transform(CRS_UTM)
regiao <- st_buffer(mrj, 40000)   # not "area": the ERA5 functions use a global `area`

# ---- 1. stations ------------------------------------------------------------------
meteo <- buscar_meteo_periodo_bulk(data_ini, data_fim, estatisticas = unique(c("median", "max", VAR)))
meteo <- meteo[!is.na(meteo$latitude) & !is.na(meteo$longitude) & !is.na(meteo[[col_est]]), ]
meteo$source <- c(ofic = "official", aura = "aurassure", wunder = "wunder")[meteo$origem]
meteo$source[meteo$estacao_id %in% AIRPORTS] <- "airport"
stations <- data.frame(station_id = meteo$estacao_id, source = meteo$source,
                       lon = meteo$longitude, lat = meteo$latitude)
stations <- stations[!duplicated(stations$station_id), ]
obs <- data.frame(station_id = meteo$estacao_id, date = as.Date(meteo$data), value = meteo[[col_est]])
obs <- obs[obs$date >= data_ini & obs$date <= data_fim, ]
sources <- data.frame(source = c("official", "airport", "aurassure", "wunder"),
                      bias = c(FALSE, TRUE, TRUE, TRUE), obs_sd = c(0.3, 0.3, NA, NA))
sources <- sources[sources$source %in% stations$source, ]
cat("\nEstacoes por rede:\n"); print(table(stations$source))

# ---- ERA5-Land, cells within 40 km of the city (as in the article) -------------------------
download_era5_days(ANO, MES, days = seq_len(as.integer(format(data_fim, "%d"))),
                   variavel = "2m_temperature", prefix = "t2m", dir_out = DIR_ERA5)
arq_nc <- list.files(DIR_ERA5, sprintf("^t2m.*%d%02d.*\\.nc$", ANO, MES), full.names = TRUE)
if (!length(arq_nc)) arq_nc <- list.files(DIR_ERA5, "^t2m_.*\\.nc$", full.names = TRUE)
era5 <- do.call(rbind, lapply(arq_nc, nc_to_df, varname = var_era5, stat = VAR))
era5 <- era5[as.Date(era5$date) >= data_ini & as.Date(era5$date) <= data_fim, ]
if (!nrow(era5)) stop("Sem ERA5 para ", tag, " em ", DIR_ERA5, call. = FALSE)
dentro <- st_intersects(st_transform(st_as_sf(era5, coords = c("lon", "lat"), crs = 4326), CRS_UTM),
                        regiao, sparse = FALSE)[, 1]
background <- data.frame(lon = era5$lon[dentro], lat = era5$lat[dentro],
                         date = as.Date(era5$date[dentro]), value = era5[[var_era5]][dentro])

# ---- 2. diagnostic: ERA5 minus station, and optional quality control --------------------
st_utm <- st_coordinates(st_transform(st_as_sf(stations, coords = c("lon", "lat"), crs = 4326), CRS_UTM))
bg_utm <- st_coordinates(st_transform(st_as_sf(background, coords = c("lon", "lat"), crs = 4326), CRS_UTM))
cel <- unique(data.frame(lon = background$lon, lat = background$lat, X = bg_utm[, 1], Y = bg_utm[, 2]))
perto <- vapply(seq_len(nrow(st_utm)), function(i) which.min((cel$X - st_utm[i, 1])^2 + (cel$Y - st_utm[i, 2])^2), 1L)
cmp <- merge(merge(obs, data.frame(station_id = stations$station_id, source = stations$source,
                                   lon = cel$lon[perto], lat = cel$lat[perto]), by = "station_id"),
             background, by = c("lon", "lat", "date"), suffixes = c("_est", "_era5"))
cmp$dif <- cmp$value_era5 - cmp$value_est
cat("\n=== ERA5 menos estacao (C), celula mais proxima, mesmo dia ===\n")
print(round(do.call(rbind, lapply(split(cmp$dif, cmp$source), function(v) c(n = length(v), media = mean(v), sd = sd(v)))), 2))
por_est <- do.call(rbind, lapply(split(cmp, cmp$station_id), function(z)
  data.frame(station_id = z$station_id[1], source = z$source[1], n = nrow(z),
             media = round(mean(z$dif), 2), sd = round(sd(z$dif), 2))))
por_est$desvio_rede <- round(por_est$media - tapply(por_est$media, por_est$source, median)[por_est$source], 2)
por_est <- por_est[order(-abs(por_est$desvio_rede)), ]
cat("\nEstacoes mais distantes da propria rede:\n"); print(head(por_est, 10), row.names = FALSE)
if (!is.na(QC)) {
  ref_src <- sources$source[!sources$bias]
  ruins <- por_est$station_id[!por_est$source %in% ref_src &
                                (por_est$n < 5 | abs(por_est$desvio_rede) > QC | (!is.na(por_est$sd) & por_est$sd > QC))]
  cat(sprintf("Controle de qualidade (qc=%g): %d fora: %s\n", QC, length(ruins),
              if (length(ruins)) paste(ruins, collapse = ", ") else "-"))
  obs <- obs[!obs$station_id %in% ruins, ]
  stations <- stations[!stations$station_id %in% ruins, ]
}
write.csv(por_est, file.path(DIR, "diagnostico_era5_estacoes.csv"), row.names = FALSE)

# ---- covariates (as in the article; cached) -----------------------------------------------
arq_cov <- file.path(DIR_COV, c(elev = "dem.tif", urb = "urb.tif", svf = "svf.tif"))
names(arq_cov) <- c("elev", "urb", "svf")
if (!all(file.exists(arq_cov))) {
  cat("\nConstruindo covariaveis (so na primeira vez)...\n")
  dem <- elevatr::get_elev_raster(locations = st_transform(regiao, 4326), z = 9, clip = "locations") |>
    rast() |> project(paste0("EPSG:", CRS_UTM))
  dem[dem < 0] <- 0
  uso <- read_sf(USO_SOLO)
  classes_urbanas <- c("Áreas residenciais", "Áreas de comércio e serviços",
                       "Áreas industriais", "Áreas de transporte",
                       "Áreas institucionais e de infraestrutura pública",
                       "Áreas de educação e saúde", "Favela")
  uso$urbano <- as.integer(uso$usoagregad %in% classes_urbanas)
  r_ref <- rast(vect(mrj), resolution = 5, crs = paste0("EPSG:", CRS_UTM))
  r_urb <- rasterize(vect(uso), r_ref, field = "urbano", fun = "mean", background = 0) |>
    aggregate(fact = 100, fun = mean, na.rm = TRUE)
  rm(r_ref, uso); gc()
  lcz_map <- LCZ4r::lcz_get_map2(LCZ4r::lcz_get_map_generator(ID = "3110e623fbe4e73b1cde55f0e9832c4f5640ac21"),
                                 city = "Rio de Janeiro")
  svf <- LCZ4r::lcz_get_parameters(lcz_map, istack = FALSE, iselect = "SVFmean")
  writeRaster(dem, arq_cov[["elev"]], overwrite = TRUE)
  writeRaster(r_urb, arq_cov[["urb"]], overwrite = TRUE)
  writeRaster(svf, arq_cov[["svf"]], overwrite = TRUE)
}
covs <- lapply(arq_cov, rast)

d <- bruce_data(obs, stations, sources, crs = CRS_UTM, background = background, background_sd = 0.7)
d <- bruce_covariates(d, covs,
                      background_method = c(elev = "buffer_mean", urb = "buffer_mean", svf = "zero"),
                      background_buffer = 5500, na_value = c(elev = 0, urb = 0))
print(d)

# ---- 3. fit ----------------------------------------------------------------------------
arq_fit <- file.path(DIR, sprintf("bruce_%s.rds", tag))
if (file.exists(arq_fit)) {
  fit <- readRDS(arq_fit); cat("\nAjuste ja existia: carregado (apague", arq_fit, "para refazer)\n")
} else {
  mesh <- bruce_mesh(d, max_edge = c(5000, 40000), cutoff = 3000, offset = c(1000, 50000))
  spec <- bruce_spec(range_prior = 15000, intercept_prior = c(0, 0.001), int_strategy = ESTRATEGIA)
  fit <- bruce_fit(d, mesh, spec, window = c(data_ini, data_fim), keep_config = TRUE)
  saveRDS(fit, arq_fit)
}
cat("\n=== Ajuste ===\n")
summary(fit)                       # includes bruce_check()

# ---- 4. maps -----------------------------------------------------------------------------
grade <- bruce_grid(mrj, RESOLUCAO, crs = fit)
rs <- c(elev = if (RESOLUCAO >= 250) "average" else "bilinear", urb = "near", svf = "near")
pred <- bruce_predict(fit, grid = grade, n_samples = N_AMOSTRAS, cap = FALSE, resample = rs, seed = NULL)
writeRaster(pred$mean, file.path(DIR, sprintf("media_%s.tif", tag)), overwrite = TRUE)
writeRaster(pred$sd, file.path(DIR, sprintf("incerteza_sd_%s.tif", tag)), overwrite = TRUE)
dias <- names(pred$mean)
media_cidade <- global(pred$mean, "mean", na.rm = TRUE)[, 1]
dia_quente <- dias[which.max(media_cidade)]

grade_paineis <- function(r, arq, cores, titulo, lim = NULL) {
  nc <- 6; nl <- ceiling(nlyr(r) / nc)
  if (is.null(lim)) lim <- quantile(values(r), c(0.005, 0.995), na.rm = TRUE)
  png(arq, width = 2400, height = 330 * nl + 60, res = 110)
  par(mfrow = c(nl, nc), mar = c(0.5, 0.5, 1.6, 3.5), oma = c(0, 0, 2, 0))
  for (i in seq_len(nlyr(r))) {
    plot(r[[i]], range = lim, col = cores, axes = FALSE, main = names(r)[i],
         legend = i %% nc == 0 | i == nlyr(r))
  }
  mtext(titulo, outer = TRUE, cex = 1.2)
  dev.off()
}
grade_paineis(pred$mean, file.path(DIR, sprintf("mapas_dias_%s.png", tag)),
              hcl.colors(50, "YlOrRd", rev = TRUE), sprintf("Temperatura (%s) por dia, mesma escala (C)", VAR))
grade_paineis(pred$sd, file.path(DIR, sprintf("incerteza_dias_%s.png", tag)),
              hcl.colors(50, "Purples", rev = TRUE), "Incerteza (sd, C)")

dec <- bruce_decompose(fit, grid = grade, dates = dia_quente, cap = FALSE, resample = rs)
png(file.path(DIR, sprintf("decomposicao_%s_%s.png", tag, dia_quente)), width = 1500, height = 520, res = 110)
plot(dec, date = dia_quente)
dev.off()

q <- t(sapply(seq_along(dias), function(i) quantile(values(pred$mean[[i]]), c(0.05, 0.5, 0.95), na.rm = TRUE)))
ref <- d$obs[d$obs$source == "official", ]
ref_dia <- tapply(ref$value, ref$date, mean)
x <- as.Date(dias)
png(file.path(DIR, sprintf("serie_%s.png", tag)), width = 1300, height = 550, res = 110)
par(mar = c(4, 4.5, 2, 1))
plot(x, q[, 2], type = "n", ylim = range(c(q, ref_dia), na.rm = TRUE), xlab = "", ylab = "Temperatura (C)",
     main = sprintf("%s: mapa (mediana e p05-p95 da cidade) e media das estacoes oficiais", tag))
polygon(c(x, rev(x)), c(q[, 1], rev(q[, 3])), col = adjustcolor("orange", 0.3), border = NA)
lines(x, q[, 2], lwd = 2, col = "darkorange3")
points(as.Date(names(ref_dia)), ref_dia, pch = 19, col = "steelblue")
legend("bottomleft", c("mapa: mediana", "mapa: p05-p95", "estacoes oficiais"), bty = "n",
       col = c("darkorange3", adjustcolor("orange", 0.5), "steelblue"), lwd = c(2, 8, NA), pch = c(NA, NA, 19))
dev.off()
cat(sprintf("\nDia mais quente (media da cidade): %s, %.1f C\n", dia_quente, max(media_cidade)))

# ---- 5. validation: leave-station-out on official stations ---------------------------------
if (RODAR_CV) {
  ofic <- fit$data$obs[fit$data$obs$source == "official", ]
  n_dias <- table(ofic$station_id)
  cand <- sort(names(n_dias)[n_dias >= 15])
  set.seed(42)
  ids <- sort(if (DOBRAS > 0 && length(cand) > DOBRAS) sample(cand, DOBRAS) else cand)
  cat(sprintf("\n=== LSO: %d estacoes oficiais deixadas de fora, uma de cada vez ===\n", length(ids)))
  cv <- bruce_cv(fit, scheme = "lso", fold_ids = ids, n_samples = N_AMOSTRAS, seed = NULL,
                 checkpoint_dir = file.path(DIR, "cv_lso"))
  print(cv)
  saveRDS(cv, file.path(DIR, sprintf("lso_%s.rds", tag)))
  p <- as.data.frame(cv$predictions)
  por <- do.call(rbind, lapply(split(p, p$station_id), function(z)
    data.frame(estacao = z$station_id[1], n = nrow(z), rmse = round(sqrt(mean((z$pred - z$obs)^2)), 2),
               vies = round(mean(z$pred - z$obs), 2),
               cobertura = round(mean(z$obs >= z$lower & z$obs <= z$upper), 2))))
  cat("\nPor estacao:\n"); print(por[order(-por$rmse), ], row.names = FALSE)
  write.csv(por, file.path(DIR, sprintf("lso_por_estacao_%s.csv", tag)), row.names = FALSE)

  png(file.path(DIR, sprintf("lso_%s.png", tag)), width = 1300, height = 600, res = 110)
  par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
  lim <- range(c(p$obs, p$pred, p$lower, p$upper), na.rm = TRUE)
  cores <- hcl.colors(length(unique(p$station_id)), "Dark 3")[as.integer(factor(p$station_id))]
  plot(p$obs, p$pred, pch = 19, cex = 0.7, col = cores, xlim = lim, ylim = lim,
       xlab = "Medido na estacao deixada de fora (C)", ylab = "Previsto pelo bruce (C)",
       main = sprintf("LSO: RMSE %.2f C, vies %+.2f C", sqrt(mean((p$pred - p$obs)^2)), mean(p$pred - p$obs)))
  segments(p$obs, p$lower, p$obs, p$upper, col = adjustcolor(cores, 0.25))
  abline(0, 1, lty = 2)
  barplot(sort(setNames(por$rmse, por$estacao)), horiz = TRUE, las = 1, cex.names = 0.7,
          xlab = "RMSE (C)", main = "Erro por estacao")
  dev.off()
}

if (RODAR_LRO) {
  ap <- st_read(MALHA_AP, quiet = TRUE)
  ap$regiao <- substr(ap$AP_RESID, 1, 1)
  cat("\n=== LRO: regioes de planejamento deixadas de fora ===\n")
  cv_lro <- bruce_cv(fit, scheme = "lro", regions = ap, region_col = "regiao", n_samples = N_AMOSTRAS,
                     seed = NULL, checkpoint_dir = file.path(DIR, "cv_lro"))
  print(cv_lro)
  saveRDS(cv_lro, file.path(DIR, sprintf("lro_%s.rds", tag)))
}

cat("\n=== Pronto. Arquivos em", DIR, "===\n")
print(list.files(DIR, pattern = "\\.(png|tif|csv|rds)$"))
