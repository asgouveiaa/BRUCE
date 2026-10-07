# ==================================================================
# 12_orvalho_ajuste.R
# First bruce fit of daily DEW POINT (Td) for Rio de Janeiro, one month.
#
# Networks (query_orvalho_periodo.R):
#   inmet    3 city stations, Td measured           -> reference, sd 0.3 C
#   redemet  airports, Td from METAR (whole degrees) -> reference, sd 0.5 C
#   ar       Alerta Rio, Td from T + RH             -> offset and noise estimated
#   aura     Aurassure, Td from T + RH              -> offset and noise estimated
#   wunder   Wunderground, Td from T + RH           -> offset and noise estimated
# Background: ERA5-Land 2 m dew point (d2m), sd 0.7 C, cells within 40 km.
# Covariates: elevation, urban fraction and distance to the sea (computed from
# the DEM: cells at or below 0 m are sea). Sky view factor is left out: it
# matters for heat retention, not for humidity.
#
# The Alerta Rio sensors look humid-biased (Td ~1.7 C above the others and
# saturated in ~40% of days in jan/2026); its offset should come out large
# and positive. MAX_FRAC_SATURADA < 1 drops station-days with more than that
# share of saturated hours (sensitivity check).
#
# Variants can be chosen on the command line and run in parallel (each with
# its own results, log and LSO checkpoints):
#   sat=0.25      drop station-days with more than 25% saturated hours
#   sigma=2       PC prior on the field sd: P(sd > 2 C) = 0.01 (default 8)
#   sem_ar=1      leave Alerta Rio out
#   range=50000 range_p=0.95   PC prior on the range: P(range < 50 km) = 0.95
#                 (default 15000 and 0.5, as in the temperature model)
#   bg_sd=2       background (ERA5) sd: trust ERA5 less (default 0.7)
#   sem_bg=1      no ERA5 background at all (stations only)
#   buffer=10     ERA5 cells within this many km of the city (default 40)
#                 The last three test whether ERA5 is what pushes the range
#                 to hundreds of km.
#   aresta=2000 cutoff=1000   finer mesh near the data (default 5000 / 3000 m)
#   dist_log=1    distance to the sea as log(1 + km) instead of km
#   ano=2023 mes=11   month to fit (default 2026-01)
#   colapso=fix   if the field collapses, refit with the day-to-day field
#                 correlation fixed (as the temperature article did);
#                 default warn = only report
#   grho=0.9      fix that correlation from the start
#   excluir=wunder:IRIODE84,wunder:IRIODEJA96   leave these stations out
#   qc=3          drop citizen stations whose difference to ERA5 is more than
#                 3 C away from their network's typical difference, or that
#                 swing more than 3 C day to day (or have < 5 days)
#   so_diagnostico=1  only print ERA5 minus station (seconds) and stop
#   threads=5:1   INLA threads (when several run at once, split the cores)
#   so_preparar=1 only query the stations, download ERA5 and build the
#                 covariates (cache), then stop: run this once before
#                 launching variants in parallel
# The station data and covariates are prepared once and cached in
# resultados/orvalho/cache, so parallel runs do not query the databases.
#
# Needs bruce >= 0.1.0.9. Rode de dentro de analises/rio:
#   Rscript 12_orvalho_ajuste.R > log_orvalho.log 2>&1
#   Rscript 12_orvalho_ajuste.R sigma=2 threads=5:1 > log_orvalho_sig2.log 2>&1 &
# ==================================================================

library(bruce)
library(terra)
library(sf)
source("config.R")
source("query_orvalho_periodo.R")
source(ERA5_FUNCOES)
if (packageVersion("bruce") < "0.1.0.9") stop("Precisa do bruce >= 0.1.0.9 (bruce_cv(fit))")

VAR <- "median"       # "median", "max" or "min" of the daily Td
args <- commandArgs(trailingOnly = TRUE)
arg <- function(nome, padrao) {
  a <- grep(paste0("^", nome, "="), args, value = TRUE)
  if (length(a)) sub("^[^=]*=", "", a[1]) else padrao
}
ANO <- as.integer(arg("ano", 2026)); MES <- as.integer(arg("mes", 1))
MAX_FRAC_SATURADA <- as.numeric(arg("sat", 1))   # 1 = keep all; e.g. 0.25 for the sensitivity check
SIGMA_U <- as.numeric(arg("sigma", 8))            # PC prior on the field sd: P(sd > SIGMA_U) = 0.01
SEM_ALERTARIO <- arg("sem_ar", "0") == "1"
RANGE_R0 <- as.numeric(arg("range", 15000))       # PC prior on the range: P(range < RANGE_R0) = RANGE_P
RANGE_P  <- as.numeric(arg("range_p", 0.5))
BG_SD <- as.numeric(arg("bg_sd", 0.7))
SEM_BG <- arg("sem_bg", "0") == "1"
BUFFER_KM <- as.numeric(arg("buffer", 40))
DIST_LOG <- arg("dist_log", "0") == "1"
EXCLUIR <- strsplit(arg("excluir", ""), ",")[[1]]
QC <- as.numeric(arg("qc", NA))
COLAPSO <- arg("colapso", "warn")
GRHO <- as.numeric(arg("grho", NA))
ARESTA <- as.numeric(arg("aresta", 5000))
CUTOFF <- as.numeric(arg("cutoff", 3000))
THREADS <- arg("threads", NA)
RODAR_LSO <- TRUE
BLOCO_DIAS_ERA5 <- 7      # if the whole-month ERA5 request fails: days per request (1 = one day at a time)
CRS_UTM <- 31983
DIR_OUT <- file.path(RES, "orvalho")
dir.create(DIR_OUT, showWarnings = FALSE, recursive = TRUE)
tag <- paste0(sprintf("td_%s_%d-%02d", VAR, ANO, MES),
              if (MAX_FRAC_SATURADA < 1) sprintf("_sat%02d", round(100 * MAX_FRAC_SATURADA)) else "",
              if (SIGMA_U != 8) sprintf("_sig%g", SIGMA_U) else "",
              if (SEM_ALERTARIO) "_semAR" else "",
              if (RANGE_R0 != 15000 || RANGE_P != 0.5) sprintf("_rng%gk_p%g", RANGE_R0 / 1000, RANGE_P) else "",
              if (BG_SD != 0.7) sprintf("_bgsd%g", BG_SD) else "",
              if (SEM_BG) "_semBG" else "",
              if (BUFFER_KM != 40) sprintf("_buf%gkm", BUFFER_KM) else "",
              if (DIST_LOG) "_logmar" else "",
              if (ARESTA != 5000 || CUTOFF != 3000) sprintf("_malha%g_%g", ARESTA, CUTOFF) else "",
              if (COLAPSO != "warn") paste0("_col", COLAPSO) else "",
              if (!is.na(GRHO)) sprintf("_grho%g", GRHO) else "",
              if (length(EXCLUIR)) sprintf("_sem%d", length(EXCLUIR)) else "",
              if (!is.na(QC)) sprintf("_qc%g", QC) else "")
cat("Variante:", tag, "\n")
DIR_CACHE <- file.path(DIR_OUT, "cache")
dir.create(DIR_CACHE, showWarnings = FALSE)
data_ini <- as.Date(sprintf("%d-%02d-01", ANO, MES))
data_fim <- seq(data_ini, by = "month", length.out = 2)[2] - 1

mrj  <- geobr::read_municipality(code_muni = 3304557, year = 2020) |> st_transform(CRS_UTM)
regiao <- st_buffer(mrj, 40000)   # not "area": the ERA5 functions use a global `area`

# ---- stations ---------------------------------------------------------------
# Queried once (no saturation filter) and cached; the variants filter here
arq_td <- file.path(DIR_CACHE, sprintf("td_obs_%s_%d-%02d.rds", VAR, ANO, MES))
if (file.exists(arq_td)) {
  td <- readRDS(arq_td)
} else {
  td <- buscar_orvalho_periodo(data_ini, data_fim, estatisticas = VAR, max_frac_saturada = 1)
  saveRDS(td, arq_td)
}
td <- td[is.na(td$frac_saturada) | td$frac_saturada <= MAX_FRAC_SATURADA, ]
if (SEM_ALERTARIO) td <- td[td$origem != "ar", ]
col <- paste0(VAR, "_orvalho")
td <- td[!is.na(td$latitude) & !is.na(td$longitude) & !is.na(td[[col]]), ]
td$source <- c(inmet = "inmet", redemet = "airport", ar = "alerta_rio", aura = "aurassure",
               wunder = "wunder")[td$origem]
# Alerta Rio and Aurassure both use numeric ids: prefix with the network so they never collide
td$estacao_id <- paste(td$origem, td$estacao_id, sep = ":")
if (length(EXCLUIR)) {
  fora <- intersect(EXCLUIR, td$estacao_id)
  cat("Excluidas:", if (length(fora)) paste(fora, collapse = ", ") else "nenhuma (ids nao encontrados)", "\n")
  td <- td[!td$estacao_id %in% EXCLUIR, ]
}

stations <- unique(data.frame(station_id = td$estacao_id, source = td$source,
                              lon = td$longitude, lat = td$latitude))
stations <- stations[!duplicated(stations$station_id), ]
obs <- data.frame(station_id = td$estacao_id, date = as.Date(td$data), value = td[[col]])
sources <- data.frame(
  source = c("inmet", "airport", "alerta_rio", "aurassure", "wunder"),
  bias   = c(FALSE,   FALSE,     TRUE,        TRUE,        TRUE),
  obs_sd = c(0.3,     0.5,       NA,          NA,          NA))
sources <- sources[sources$source %in% stations$source, ]
cat("\nEstacoes por rede:\n"); print(table(stations$source))

# ---- background: ERA5-Land dew point -----------------------------------------
# Whole month first; if the Copernicus request fails, weekly blocks (each its
# own file, joined afterwards). Files already downloaded are not requested again.
dias_mes <- seq_len(as.integer(format(data_fim, "%d")))
var_era5 <- paste0("d2m_", VAR)
arq_nc <- file.path(DIR_ERA5, sprintf("d2m_%d%02d.nc", ANO, MES))
download_era5_days(ANO, MES, days = dias_mes, variavel = "2m_dewpoint_temperature",
                   prefix = "d2m", dir_out = DIR_ERA5)
if (file.exists(arq_nc)) {
  era5 <- nc_to_df(arq_nc, varname = var_era5, stat = VAR)
} else {
  message(sprintf("Pedido do mes inteiro falhou; baixando o d2m em blocos de %d dia(s)", BLOCO_DIAS_ERA5))
  blocos <- split(dias_mes, ceiling(dias_mes / BLOCO_DIAS_ERA5))
  era5 <- do.call(rbind, lapply(seq_along(blocos), function(b) {
    pref <- sprintf("d2m_%dd_b%d", BLOCO_DIAS_ERA5, b)
    download_era5_days(ANO, MES, days = blocos[[b]], variavel = "2m_dewpoint_temperature",
                       prefix = pref, dir_out = DIR_ERA5)
    f <- file.path(DIR_ERA5, sprintf("%s_%d%02d.nc", pref, ANO, MES))
    if (!file.exists(f)) {
      stop("O ERA5 d2m de ", ANO, "-", sprintf("%02d", MES), " (dias ", min(blocos[[b]]), "-",
           max(blocos[[b]]), ") nao foi baixado. Veja o motivo em ",
           "https://cds.climate.copernicus.eu/requests e rode de novo; o que ja foi baixado e pulado.",
           call. = FALSE)
    }
    nc_to_df(f, varname = var_era5, stat = VAR)
  }))
}
inside <- st_intersects(st_transform(st_as_sf(era5, coords = c("lon", "lat"), crs = 4326), CRS_UTM),
                        st_buffer(mrj, BUFFER_KM * 1000), sparse = FALSE)[, 1]
background <- data.frame(lon = era5$lon[inside], lat = era5$lat[inside],
                         date = as.Date(era5$date[inside]), value = era5[[var_era5]][inside])

# ---- covariates ------------------------------------------------------------------
# Built once and cached as GeoTIFF (they do not depend on the month)
arq_cov <- file.path(DIR_CACHE, c(dem = "dem.tif", dist_mar = "dist_mar.tif", urb = "urb.tif"))
names(arq_cov) <- c("dem", "dist_mar", "urb")
if (!all(file.exists(arq_cov))) {
  dem_bruto <- elevatr::get_elev_raster(locations = st_transform(regiao, 4326), z = 9, clip = "locations") |>
    rast() |> project(paste0("EPSG:", CRS_UTM))
  mar <- ifel(dem_bruto <= 0, 1, NA)              # sea (and bay) cells
  dist_mar <- distance(mar) / 1000                 # km to the nearest sea cell
  names(dist_mar) <- "dist_mar"
  dem <- dem_bruto; dem[dem < 0] <- 0

  uso <- read_sf(USO_SOLO)
  classes_urbanas <- c("Áreas residenciais", "Áreas de comércio e serviços", "Áreas industriais",
                       "Áreas de transporte", "Áreas institucionais e de infraestrutura pública",
                       "Áreas de educação e saúde", "Favela")
  uso$urbano <- as.integer(uso$usoagregad %in% classes_urbanas)
  r_ref <- rast(vect(mrj), resolution = 5, crs = paste0("EPSG:", CRS_UTM))
  r_urb <- rasterize(vect(uso), r_ref, field = "urbano", fun = "mean", background = 0) |>
    aggregate(fact = 100, fun = mean, na.rm = TRUE)
  rm(r_ref, uso); gc()
  writeRaster(dem, arq_cov[["dem"]], overwrite = TRUE)
  writeRaster(dist_mar, arq_cov[["dist_mar"]], overwrite = TRUE)
  writeRaster(r_urb, arq_cov[["urb"]], overwrite = TRUE)
}
dem <- rast(arq_cov[["dem"]]); dist_mar <- rast(arq_cov[["dist_mar"]]); r_urb <- rast(arq_cov[["urb"]])
if (DIST_LOG) dist_mar <- log1p(dist_mar)   # first kilometres count more
if (arg("so_preparar", "0") == "1") {
  cat("Cache pronto em", DIR_CACHE, "- agora as variantes podem rodar em paralelo\n")
  quit(save = "no")
}

# ---- diagnostic: ERA5 against the reference stations (nearest cell, same day) --
st_utm <- st_coordinates(st_transform(st_as_sf(stations, coords = c("lon", "lat"), crs = 4326), CRS_UTM))
bg_utm <- st_coordinates(st_transform(st_as_sf(background, coords = c("lon", "lat"), crs = 4326), CRS_UTM))
cel <- unique(data.frame(lon = background$lon, lat = background$lat, X = bg_utm[, 1], Y = bg_utm[, 2]))
mais_perto <- vapply(seq_len(nrow(st_utm)), function(i)
  which.min((cel$X - st_utm[i, 1])^2 + (cel$Y - st_utm[i, 2])^2), 1L)
st_cel <- data.frame(station_id = stations$station_id, source = stations$source,
                     lon = cel$lon[mais_perto], lat = cel$lat[mais_perto])
cmp <- merge(merge(obs, st_cel, by = "station_id"), background, by = c("lon", "lat", "date"),
             suffixes = c("_est", "_era5"))
cat("\n=== ERA5 menos estacao (C), celula mais proxima, mesmo dia ===\n")
print(round(do.call(rbind, lapply(split(cmp, cmp$source), function(z)
  c(n = nrow(z), media = mean(z$value_era5 - z$value_est), sd = sd(z$value_era5 - z$value_est)))), 2))
cat("Celulas ERA5 por dia:", nrow(cel), "  estacoes:", nrow(stations), "\n")
por_est <- do.call(rbind, lapply(split(cmp, cmp$station_id), function(z)
  data.frame(station_id = z$station_id[1], source = z$source[1], n = nrow(z),
             media = round(mean(z$value_era5 - z$value_est), 2), sd = round(sd(z$value_era5 - z$value_est), 2))))
# distance from the TYPICAL difference of the station's own network (the
# network's bias is expected; a station far from its own network is suspect)
med_rede <- tapply(por_est$media, por_est$source, median)
por_est$desvio_rede <- round(por_est$media - med_rede[por_est$source], 2)
por_est <- por_est[order(-abs(por_est$desvio_rede)), ]
cat("\nEstacoes mais distantes da propria rede (ERA5 - estacao, menos a mediana da rede):\n")
print(head(por_est, 10), row.names = FALSE)
# automatic quality control: only biased (citizen) networks, never the references
if (!is.na(QC)) {
  ref_src <- sources$source[!sources$bias]
  ruins <- por_est$station_id[!por_est$source %in% ref_src & por_est$n >= 5 &
                                (abs(por_est$desvio_rede) > QC | por_est$sd > QC)]
  ruins <- union(ruins, por_est$station_id[!por_est$source %in% ref_src & por_est$n < 5])
  cat(sprintf("\nControle de qualidade (qc=%g): %d estacao(oes) fora: %s\n", QC, length(ruins),
              if (length(ruins)) paste(ruins, collapse = ", ") else "-"))
  obs <- obs[!obs$station_id %in% ruins, ]
  stations <- stations[!stations$station_id %in% ruins, ]
}
if (arg("so_diagnostico", "0") == "1") quit(save = "no")

d <- bruce_data(obs, stations, sources, crs = CRS_UTM,
                background = if (SEM_BG) NULL else background, background_sd = BG_SD)
d <- bruce_covariates(d, list(elev = dem, urb = r_urb, dist_mar = dist_mar),
                      background_method = "buffer_mean", background_buffer = 5500,
                      na_value = c(elev = 0, urb = 0))
d

# ---- fit -------------------------------------------------------------------------
mesh <- bruce_mesh(d, max_edge = c(ARESTA, 40000), cutoff = CUTOFF, offset = c(1000, 50000))
spec <- bruce_spec(range_prior = RANGE_R0, range_prior_prob = RANGE_P,
                   intercept_prior = c(0, 0.001), int_strategy = "eb",
                   sigma_prior = c(SIGMA_U, 0.01), on_collapse = COLAPSO,
                   group_rho_fixed = if (is.na(GRHO)) NULL else GRHO,
                   num_threads = if (is.na(THREADS)) NULL else THREADS)
fit <- bruce_fit(d, mesh, spec, window = c(data_ini, data_fim), keep_config = TRUE)
summary(fit)

cat("\n=== Vies de cada rede em relacao a referencia (C) ===\n")
fx <- fit$inla$summary.fixed
print(round(fx[grep("^bias_", rownames(fx)), c("mean", "0.025quant", "0.975quant")], 2))

grade <- bruce_grid(mrj, 500, crs = fit)
rs <- c(elev = "average", urb = "near", dist_mar = "average")
pred <- bruce_predict(fit, grid = grade, resample = rs)
dec  <- bruce_decompose(fit, grid = grade, resample = rs)
print(dec)
dia <- names(which.max(tapply(obs$value, obs$date, mean)))
png(file.path(DIR_OUT, sprintf("decomposicao_%s.png", tag)), width = 1500, height = 520, res = 110)
plot(dec, date = dia)
dev.off()

saveRDS(list(fit = fit, pred_mean = wrap(pred$mean), variance = dec$variance),
        file.path(DIR_OUT, sprintf("bruce_%s.rds", tag)))
writeRaster(pred$mean, file.path(DIR_OUT, sprintf("bruce_pred_%s.tif", tag)), overwrite = TRUE)
message("Ajuste salvo em ", DIR_OUT)

# ---- leave-station-out on the reference stations (INMET + airports) ----------------
if (RODAR_LSO) {
  cv <- bruce_cv(fit, n_samples = 100, seed = NULL,
                 checkpoint_dir = file.path(DIR_OUT, sprintf("cv_lso_%s", tag)))
  print(cv)
  saveRDS(cv, file.path(DIR_OUT, sprintf("bruce_lso_%s.rds", tag)))
}
