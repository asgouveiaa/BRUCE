# ==================================================================
# rio_adapter.R
# Maps the Rio de Janeiro inputs of modelo_9meses_median_max.R onto the
# bruce data contract and fits one month. Everything Rio-specific lives
# here; the package itself knows nothing about INMET, Aurassure or ERA5.
#
# Correspondence with the original script:
#   origem == "ofic"           -> source "official"  (bias = FALSE, obs_sd = 0.3)
#   ESTACOES_AEROPORTO         -> source "airport"   (bias = TRUE,  obs_sd = 0.3)
#   origem == "aura"           -> source "aurassure" (bias = TRUE,  obs_sd = NA)
#   origem == "wunder"         -> source "wunder"    (bias = TRUE,  obs_sd = NA)
#   ERA5-Land t2m              -> background         (sd = 0.7)
#   elev, frac_urb (5.5 km buffer mean at ERA5 cells), svf (0 at ERA5 cells)
#
# The airport offset used to share the official likelihood; here airports
# get their own likelihood with the same fixed sd, which is the same model.
# ==================================================================

library(bruce)
library(terra)
library(sf)

source("~/scratch/spde_paper/query_meteorj_periodo.R")
source("~/MEGA/cisc/era5_interp/com_aura_novo/funcoes downscaling ERA5 MRJ.R")

ANO <- 2025; MES <- 1; VAR <- "median"     # "median", "max" or "min"
CRS_UTM <- 31983
AIRPORTS <- c("SBJR", "SBGL", "SBAF", "SBSC")

col_est  <- paste0(VAR, "_temperatura")   # median_temperatura / max_temperatura / min_temperatura
var_era5 <- paste0("t2m_", VAR)           # t2m_median / t2m_max / t2m_min

data_ini <- as.Date(sprintf("%d-%02d-01", ANO, MES))
data_fim <- seq(data_ini, by = "month", length.out = 2)[2] - 1

# ---- stations -------------------------------------------------------------
meteo <- buscar_meteo_periodo_bulk(data_ini, data_fim)
meteo <- meteo[!is.na(meteo$latitude) & !is.na(meteo[[col_est]]), ]
meteo$source <- c(ofic = "official", aura = "aurassure", wunder = "wunder")[meteo$origem]
meteo$source[meteo$estacao_id %in% AIRPORTS] <- "airport"

stations <- unique(data.frame(station_id = meteo$estacao_id, source = meteo$source,
                              lon = meteo$longitude, lat = meteo$latitude))
stations <- stations[!duplicated(stations$station_id), ]
obs <- data.frame(station_id = meteo$estacao_id, date = meteo$data, value = meteo[[col_est]])

sources <- data.frame(
  source = c("official", "airport", "aurassure", "wunder"),
  bias   = c(FALSE,      TRUE,      TRUE,        TRUE),
  obs_sd = c(0.3,        0.3,       NA,          NA)
)

# ---- background (ERA5-Land) -----------------------------------------------
DIR_ERA5 <- "~/MEGA/cisc/era5_interp/independente"
download_era5_days(ANO, MES, days = seq_len(as.integer(format(data_fim, "%d"))),
                   variavel = "2m_temperature", prefix = "t2m", dir_out = DIR_ERA5)
era5 <- do.call(rbind, lapply(list.files(DIR_ERA5, "^t2m_.*\\.nc$", full.names = TRUE),
                              nc_to_df, varname = var_era5, stat = VAR))
era5 <- era5[era5$date >= data_ini & era5$date <= data_fim, ]
mrj <- geobr::read_municipality(code_muni = 3304557, year = 2020) |> st_transform(CRS_UTM)
area <- st_buffer(mrj, 40000)
keep <- st_intersects(st_transform(st_as_sf(era5, coords = c("lon", "lat"), crs = 4326), CRS_UTM),
                      area, sparse = FALSE)[, 1]
background <- data.frame(lon = era5$lon[keep], lat = era5$lat[keep],
                         date = era5$date[keep], value = era5[[var_era5]][keep])

# ---- covariates: one SpatRaster on a common fine grid ------------------------
# Use cached rasters (dem_utm.tif, frac_urb.tif, svf) resampled to a common
# 100 m grid; the prediction grid is chosen separately in bruce_predict().
# PATHS BELOW ARE PLACEHOLDERS -- point them at the real cache files
# (svf can be written once from lcz_get_parameters(..., iselect = "SVFmean")).
dem <- rast("~/Documents/lcz_painel/cache/dem_utm.tif"); dem[dem < 0] <- 0
urb <- rast("~/Documents/lcz_painel/cache/frac_urb.tif")
svf <- rast("~/Documents/lcz_painel/cache/svf_mean.tif")       # from LCZ4r
base <- rast(ext(vect(area)), resolution = 100, crs = paste0("EPSG:", CRS_UTM))
cov <- c(resample(dem, base, "bilinear"), resample(urb, base, "average"),
         resample(project(svf, base), base, "near"))
names(cov) <- c("elev", "urb", "svf")

# ---- bruce ----------------------------------------------------------------
d <- bruce_data(obs, stations, sources, crs = CRS_UTM, background = background,
                background_sd = 0.7)
d <- bruce_covariates(d, cov, background_method = c(elev = "buffer_mean", urb = "buffer_mean",
                                                    svf = "zero"),
                      background_buffer = 5500)
mesh <- bruce_mesh(d, max_edge = c(5000, 40000), cutoff = 3000, offset = c(1000, 50000))

spec <- bruce_spec(range_prior = 15000)   # value in the current scripts; "nn" = network-anchored
fit <- bruce_fit(d, mesh, spec, window = c(data_ini, data_fim), keep_config = TRUE)
summary(fit)

grid500 <- mask(rast(vect(mrj), resolution = 500, crs = paste0("EPSG:", CRS_UTM)) |> setValues(1),
                vect(mrj))
pred <- bruce_predict(fit, grid = grid500, n_samples = 200, resample = c(svf = "near"))

cv <- bruce_cv(d, mesh, spec, window = c(data_ini, data_fim), scheme = "lro",
               regions = st_read("./malhas/ap_malha.gpkg") |>
                 transform(regiao = substr(AP_RESID, 1, 1)),
               region_col = "regiao", checkpoint_dir = "cv_lro_2025-01")
cv
