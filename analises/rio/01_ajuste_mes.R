# ==================================================================
# 01_ajuste_mes.R
# Fits one month of the Rio de Janeiro model with bruce, reproducing
# modelo_9meses_median_max.R (the model of the SPDE article) as closely
# as the package allows. Rode de dentro de analises/rio (caminhos em config.R).
#
# Correspondence with the original script:
#   origem == "ofic" (except airports) -> source "official"  (bias = FALSE, sd 0.3)
#   ESTACOES_AEROPORTO                 -> source "airport"   (bias = TRUE,  sd 0.3)
#   origem == "aura"                   -> source "aurassure" (bias = TRUE,  sd ~ PC prior)
#   origem == "wunder"                 -> source "wunder"    (bias = TRUE,  sd ~ PC prior)
#   ERA5-Land t2m (40 km buffer)       -> background         (sd 0.7)
#   b0 -> intercept, elev_std -> elev, urb_std -> urb, svf_std -> svf,
#   aura -> bias_aurassure, wunder -> bias_wunder, aeroporto -> bias_airport,
#   id_estacao_idx -> station_idx, dia_idx -> day_idx
#
# Known small differences from the original (expected to move results by
# much less than one posterior sd):
#   * urban fraction at ERA5 cells is the 5.5 km buffer mean everywhere
#     (the original set it to 0 for cells more than 1 km outside the city);
#   * the SVF standardisation uses distinct stations, not station-days.
# Covariates are read from their native rasters (bruce >= 0.1.1), as in the
# article; with inputs taken from the article's .rds the fit is identical
# (see 04_comparar_entradas_identicas.R).
#
# Runtime: about the same as one month of the original (~25-30 min).
# Run in tmux:  Rscript 01_ajuste_mes.R > log_bruce_2025-01.log 2>&1
# ==================================================================

library(bruce)
library(terra)
library(sf)
source("config.R")

source("query_meteorj_periodo.R")
source(ERA5_FUNCOES)

ANO <- 2025; MES <- 1; VAR <- "median"        # "median", "max" or "min"
CRS_UTM  <- 31983
AIRPORTS <- c("SBJR", "SBGL", "SBAF", "SBSC")
DIR_OUT  <- RES
tag <- sprintf("%s_%d-%02d", VAR, ANO, MES)

col_est  <- paste0(VAR, "_temperatura")
var_era5 <- paste0("t2m_", VAR)
data_ini <- as.Date(sprintf("%d-%02d-01", ANO, MES))
data_fim <- seq(data_ini, by = "month", length.out = 2)[2] - 1

mrj  <- geobr::read_municipality(code_muni = 3304557, year = 2020) |> st_transform(CRS_UTM)
area <- st_buffer(mrj, 40000)

# ---- stations -------------------------------------------------------------
# query_meteorj_periodo.R v3: "min" adds min_temperatura (median/max unchanged)
meteo <- buscar_meteo_periodo_bulk(data_ini, data_fim, estatisticas = unique(c("median", "max", VAR)))
meteo <- meteo[!is.na(meteo$latitude) & !is.na(meteo$longitude) & !is.na(meteo[[col_est]]), ]
meteo$source <- c(ofic = "official", aura = "aurassure", wunder = "wunder")[meteo$origem]
meteo$source[meteo$estacao_id %in% AIRPORTS] <- "airport"

stations <- data.frame(station_id = meteo$estacao_id, source = meteo$source,
                       lon = meteo$longitude, lat = meteo$latitude)
stations <- stations[!duplicated(stations$station_id), ]
obs <- data.frame(station_id = meteo$estacao_id, date = as.Date(meteo$data),
                  value = meteo[[col_est]])
obs <- obs[obs$date >= data_ini & obs$date <= data_fim, ]

sources <- data.frame(
  source = c("official", "airport", "aurassure", "wunder"),
  bias   = c(FALSE,      TRUE,      TRUE,        TRUE),
  obs_sd = c(0.3,        0.3,       NA,          NA)
)

# ---- background: ERA5-Land, cells within 40 km of the city -----------------
download_era5_days(ANO, MES, days = seq_len(as.integer(format(data_fim, "%d"))),
                   variavel = "2m_temperature", prefix = "t2m", dir_out = DIR_ERA5)
era5 <- do.call(rbind, lapply(list.files(DIR_ERA5, "^t2m_.*\\.nc$", full.names = TRUE),
                              nc_to_df, varname = var_era5, stat = VAR))
era5 <- era5[era5$date >= data_ini & era5$date <= data_fim, ]
inside <- st_intersects(st_transform(st_as_sf(era5, coords = c("lon", "lat"), crs = 4326), CRS_UTM),
                        area, sparse = FALSE)[, 1]
background <- data.frame(lon = era5$lon[inside], lat = era5$lat[inside],
                         date = as.Date(era5$date[inside]), value = era5[[var_era5]][inside])

# ---- covariates: built exactly as in modelo_9meses_median_max.R ------------
dem <- elevatr::get_elev_raster(locations = st_transform(area, 4326), z = 9, clip = "locations") |>
  rast() |> project(paste0("EPSG:", CRS_UTM))
dem[dem < 0] <- 0

uso <- read_sf(USO_SOLO)
classes_urbanas <- c("Áreas residenciais", "Áreas de comércio e serviços", "Áreas industriais",
                     "Áreas de transporte", "Áreas institucionais e de infraestrutura pública",
                     "Áreas de educação e saúde", "Favela")
uso$urbano <- as.integer(uso$usoagregad %in% classes_urbanas)
r_ref <- rast(vect(mrj), resolution = 5, crs = paste0("EPSG:", CRS_UTM))
r_urb <- rasterize(vect(uso), r_ref, field = "urbano", fun = "mean", background = 0) |>
  aggregate(fact = 100, fun = mean, na.rm = TRUE)                      # 500 m, as in the article
rm(r_ref); gc()

lcz_map <- LCZ4r::lcz_get_map2(
  LCZ4r::lcz_get_map_generator(ID = "3110e623fbe4e73b1cde55f0e9832c4f5640ac21"),
  city = "Rio de Janeiro")
svf <- LCZ4r::lcz_get_parameters(lcz_map, istack = FALSE, iselect = "SVFmean")

# Each covariate stays in its native raster (as in the article): the DEM in
# UTM at its own resolution, urban fraction at 500 m over the city only, and
# SVF in lon/lat from LCZ4r. Stations read each raster at their point; ERA5
# cells take the 5.5 km buffer mean (elev, urb) or nothing (svf).
# na_value = 0 reproduces the article's replace_na(..., 0) for elevation and
# urban fraction (outside the city the urban fraction is 0).
d <- bruce_data(obs, stations, sources, crs = CRS_UTM, background = background,
                background_sd = 0.7)
d <- bruce_covariates(d, list(elev = dem, urb = r_urb, svf = svf),
                      background_method = c(elev = "buffer_mean", urb = "buffer_mean", svf = "zero"),
                      background_buffer = 5500,
                      na_value = c(elev = 0, urb = 0))
d

mesh <- bruce_mesh(d, max_edge = c(5000, 40000), cutoff = 3000, offset = c(1000, 50000))

# Priors of the article. on_collapse = "warn" so the fit is never silently
# altered: if it collapses we want to know, not to have it fixed.
spec <- bruce_spec(range_prior = 15000, intercept_prior = c(0, 0.001), on_collapse = "warn")

t0 <- Sys.time()
fit <- bruce_fit(d, mesh, spec, window = c(data_ini, data_fim), keep_config = TRUE)
message(sprintf("Fit: %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
summary(fit)

# Prediction on the 500 m grid of the article, without covariate capping
# (the article's predictions were not capped).
grid500 <- rast(vect(mrj), resolution = 500, crs = paste0("EPSG:", CRS_UTM))
values(grid500) <- 1
grid500 <- mask(grid500, vect(mrj))
pred <- bruce_predict(fit, grid = grid500, cap = FALSE,
                      resample = c(elev = "average", urb = "near", svf = "near"))
dec  <- bruce_decompose(fit, grid = grid500, cap = FALSE,
                        resample = c(elev = "average", urb = "near", svf = "near"))

saveRDS(list(fit = fit, pred_mean = wrap(pred$mean), variance = dec$variance),
        file.path(DIR_OUT, sprintf("bruce_%s.rds", tag)))
writeRaster(pred$mean, file.path(DIR_OUT, sprintf("bruce_pred_%s.tif", tag)), overwrite = TRUE)
message("Saved to ", DIR_OUT, "/. Next: 02_comparar_com_artigo.R")

# ---- optional: leave-region-out with CRPS (hours: 5 regions x ~25 min) -----
# Mirrors calcular_crps_lro_2025.R: all networks removed in the held-out
# region, official stations and airports scored, 200 samples, seed 42.
RODAR_LRO <- FALSE
if (RODAR_LRO) {
  ap <- st_read(MALHA_AP, quiet = TRUE)
  ap$regiao <- substr(ap$AP_RESID, 1, 1)
  cv <- bruce_cv(d, mesh, spec, window = c(data_ini, data_fim), scheme = "lro",
                 regions = ap, region_col = "regiao",
                 test_sources = c("official", "airport"),
                 n_samples = 200, seed = 42,
                 checkpoint_dir = file.path(DIR_OUT, sprintf("cv_lro_%s", tag)))
  print(cv)
  saveRDS(cv, file.path(DIR_OUT, sprintf("bruce_lro_%s.rds", tag)))
}
