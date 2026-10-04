# ==================================================================
# berlin_adapter.R
# bruce outside Rio: daily air temperature in Berlin from
#   * DWD stations (German Weather Service)      -> reference network
#   * FAIRUrbTemp street-level sensors (Berlin)  -> biased network
#   * ERA5-Land hourly 2 m temperature           -> background
# with covariates available for any city:
#   * elevation (elevatr), built-up fraction (ESA WorldCover 10 m),
#     sky view factor (LCZ4r, global LCZ map).
#
# Daily statistics are computed the same way for every source, as in the
# Rio pipeline: hourly values first, then the daily median / max / min of
# those hourly values, on LOCAL days (Europe/Berlin).
#
# BEFORE RUNNING
#   1. FAIRUrbTemp: download the Berlin folder from
#      https://doi.org/10.48620/93247 (BORIS portal, CC-BY), unzip it and set
#      DIR_FAIR below. Files are SEF (Station Exchange Format) .tsv files; the
#      reader below finds the hourly-mean files by their content, so the
#      folder layout does not matter.
#   2. ERA5-Land: a CDS account with the key set for ecmwfr
#      (ecmwfr::wf_set_key(...)), as for the Rio ERA5 downloads.
#   3. Packages: bruce, terra, sf, ecmwfr, elevatr, LCZ4r, geodata.
#
# Rode de dentro de analises/berlin, em tmux:
#   cd ~/bruce/analises/berlin && Rscript 01_ajuste_mes.R > log_berlim.log 2>&1
# ==================================================================

library(bruce)
library(terra)
library(sf)

# FAIRUrbTemp Berlin (TU Berlin, 11 stations): Jan 2020 - Jan 2024.
# July 2022 includes the 19-20 July heatwave.
ANO <- 2022; MES <- 7; VAR <- "median"       # "median", "max" or "min"
DIR_FAIR <- "~/Downloads/SEF_output"          # unzipped Berlin folder (all subfolders are scanned)
DIR_OUT  <- "resultados"                     # relativo a analises/berlin; fora do git
CRS_UTM  <- 25833                             # ETRS89 / UTM 33N
TZ       <- "Europe/Berlin"
dir.create(DIR_OUT, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(DIR_OUT, "cache"), showWarnings = FALSE)
tag <- sprintf("%s_%d-%02d", VAR, ANO, MES)

data_ini <- as.Date(sprintf("%d-%02d-01", ANO, MES))
data_fim <- seq(data_ini, by = "month", length.out = 2)[2] - 1

# ==================================================================
# Helpers
# ==================================================================
`%||%` <- function(a, b) if (length(a)) a else b

# Hourly values -> one row per station and local day with median/max/min.
# Days with fewer than `min_horas` hourly values are dropped.
diario <- function(h, min_horas = 18) {
  h$data <- as.Date(format(h$hora, tz = TZ, usetz = FALSE))
  h <- h[!is.na(h$valor), ]
  agg <- function(f) tapply(h$valor, list(h$station_id, h$data), f)
  n <- agg(length)
  out <- data.frame(
    station_id = rep(rownames(n), times = ncol(n)),
    date = as.Date(rep(colnames(n), each = nrow(n))),
    n = as.vector(n), median = as.vector(agg(median)),
    max = as.vector(agg(max)), min = as.vector(agg(min)))
  out[!is.na(out$n) & out$n >= min_horas, ]
}

# ---- DWD hourly air temperature (open data) -----------------------------------
DWD_BASE <- "https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/hourly/air_temperature"

dwd_estacoes <- function() {
  f <- file.path(DIR_OUT, "cache", "TU_Stundenwerte_Beschreibung_Stationen.txt")
  if (!file.exists(f)) download.file(file.path(DWD_BASE, "recent", basename(f)), f, quiet = TRUE)
  ln <- readLines(f, encoding = "latin1")[-(1:2)]
  rx <- "^(\\d{5}) (\\d{8}) (\\d{8}) +(-?\\d+) +([0-9.]+) +([0-9.]+) +(.+?) {2,}(\\S.*?) +(\\S+) *$"
  m <- regmatches(ln, regexec(rx, ln))
  m <- do.call(rbind, m[lengths(m) == 10])
  data.frame(station_id = m[, 2], inicio = as.Date(m[, 3], "%Y%m%d"),
             fim = as.Date(m[, 4], "%Y%m%d"), alt = as.numeric(m[, 5]),
             lat = as.numeric(m[, 6]), lon = as.numeric(m[, 7]),
             nome = trimws(m[, 8]), stringsAsFactors = FALSE)
}

# Reads a DWD 'produkt_tu_stunde_*.txt' file -> station_id, hora (UTC), valor
ler_dwd_produkt <- function(txt) {
  x <- read.table(txt, sep = ";", header = TRUE, strip.white = TRUE,
                  colClasses = c(STATIONS_ID = "character", MESS_DATUM = "character"))
  x <- x[x$TT_TU > -999, ]                     # -999 = missing
  data.frame(station_id = sprintf("%05d", as.integer(x$STATIONS_ID)),
             hora = as.POSIXct(x$MESS_DATUM, format = "%Y%m%d%H", tz = "UTC"),
             valor = as.numeric(x$TT_TU))
}

dwd_horario <- function(id) {
  achar <- function(pasta) {
    idx <- tryCatch(readLines(file.path(DWD_BASE, pasta), warn = FALSE), error = function(e) character())
    z <- regmatches(idx, regexpr(sprintf("stundenwerte_TU_%s_[^\"]*\\.zip", id), idx))
    if (length(z)) file.path(DWD_BASE, pasta, unique(z)) else character()
  }
  urls <- c(achar("historical"), achar("recent"))
  do.call(rbind, lapply(urls, function(u) {
    zf <- file.path(DIR_OUT, "cache", basename(u))
    if (!file.exists(zf)) download.file(u, zf, quiet = TRUE, mode = "wb")
    pasta <- file.path(DIR_OUT, "cache", sub("\\.zip$", "", basename(u)))
    unzip(zf, exdir = pasta)
    ler_dwd_produkt(list.files(pasta, "^produkt_tu_stunde", full.names = TRUE)[1])
  }))
}

# ---- FAIRUrbTemp (SEF files) ---------------------------------------------------
ler_sef <- function(f) {
  ln <- readLines(f, warn = FALSE, encoding = "UTF-8")
  i_dados <- grep("^Year\\t", ln)[1]
  if (is.na(i_dados)) return(NULL)
  cab <- strsplit(ln[seq_len(i_dados - 1)], "\t")
  meta <- stats::setNames(vapply(cab, function(z) if (length(z) > 1) z[2] else NA_character_, ""),
                          vapply(cab, `[`, "", 1))
  x <- read.delim(f, skip = i_dados - 1, header = TRUE, colClasses = "character",
                  na.strings = c("NA", ""))
  hora <- as.POSIXct(sprintf("%s-%s-%s %s:%s", x$Year, x$Month, x$Day,
                             ifelse(is.na(x$Hour), "0", x$Hour), ifelse(is.na(x$Minute), "0", x$Minute)),
                     format = "%Y-%m-%d %H:%M", tz = "UTC")
  list(meta = meta,
       dados = data.frame(station_id = meta[["ID"]], hora = hora,
                          valor = suppressWarnings(as.numeric(x$Value)),
                          flag = !is.na(x$Meta) & grepl("qc", x$Meta, ignore.case = TRUE)))
}

# All hourly-mean air-temperature series under a folder: chosen by content
# (statistic 'mean' and a typical step of one hour), not by file names.
fair_horario <- function(pasta) {
  fs <- list.files(path.expand(pasta), pattern = "\\.tsv$", recursive = TRUE, full.names = TRUE)
  if (!length(fs)) stop("No .tsv files under ", pasta)
  series <- list(); vistos <- list()
  for (f in fs) {
    cab <- readLines(f, n = 12, warn = FALSE)
    stat <- sub("^Stat\\t", "", grep("^Stat\\t", cab, value = TRUE))
    vbl  <- sub("^Vbl\\t", "", grep("^Vbl\\t", cab, value = TRUE))
    vistos[[length(vistos) + 1]] <- data.frame(Vbl = vbl[1] %||% NA, Stat = stat[1] %||% NA)
    if (!length(stat) || !identical(tolower(stat[1]), "mean")) next
    s <- ler_sef(f)
    if (is.null(s) || nrow(s$dados) < 2) next
    passo <- stats::median(as.numeric(diff(sort(s$dados$hora)), units = "mins"), na.rm = TRUE)
    if (!is.finite(passo) || passo < 55 || passo > 65) next
    s$dados$lat <- as.numeric(s$meta[["Lat"]]); s$dados$lon <- as.numeric(s$meta[["Lon"]])
    series[[length(series) + 1]] <- s$dados
  }
  if (!length(series)) {
    print(table(do.call(rbind, vistos)))
    stop("No hourly-mean series found; send the table above (Vbl x Stat found) to adjust the reader")
  }
  out <- do.call(rbind, series)
  message(sprintf("FAIRUrbTemp: %d hourly series, %d flagged values removed",
                  length(series), sum(out$flag)))
  out[!out$flag, ]
}
# ==================================================================
# 1. Data
# ==================================================================
berlim <- geodata::gadm("DEU", level = 1, path = file.path(DIR_OUT, "cache"))
berlim <- st_as_sf(berlim)
berlim <- st_transform(berlim[berlim$NAME_1 == "Berlin", ], CRS_UTM)
area   <- st_buffer(berlim, 40000)

# ---- FAIRUrbTemp ----------------------------------------------------------------
fair_h <- fair_horario(DIR_FAIR)
cobertura <- tapply(fair_h$station_id, format(fair_h$hora, "%Y-%m"), function(z) length(unique(z)))
cat("\nFAIRUrbTemp Berlin -- estacoes com dados por mes:\n"); print(cobertura)
fair_h <- fair_h[fair_h$hora >= as.POSIXct(data_ini - 1, tz = "UTC") &
                 fair_h$hora <  as.POSIXct(data_fim + 2, tz = "UTC"), ]
if (length(unique(fair_h$station_id)) < 5) {
  stop(sprintf("Only %d FAIRUrbTemp stations in %d-%02d; pick a month from the table above",
               length(unique(fair_h$station_id)), ANO, MES))
}
fair_st <- unique(fair_h[, c("station_id", "lat", "lon")])
fair_st <- fair_st[!duplicated(fair_st$station_id), ]
fair_d  <- diario(fair_h)

# ---- DWD: stations within the 40 km buffer, active in the month ----------------
dwd_all <- dwd_estacoes()
dwd_pts <- st_transform(st_as_sf(dwd_all, coords = c("lon", "lat"), crs = 4326, remove = FALSE), CRS_UTM)
dwd_sel <- dwd_all[st_intersects(dwd_pts, area, sparse = FALSE)[, 1] &
                   dwd_all$inicio <= data_ini & dwd_all$fim >= data_fim, ]
cat(sprintf("\nDWD: %d estacoes no buffer de 40 km ativas no mes: %s\n",
            nrow(dwd_sel), paste(dwd_sel$nome, collapse = ", ")))
dwd_h <- do.call(rbind, lapply(dwd_sel$station_id, dwd_horario))
dwd_h <- dwd_h[dwd_h$hora >= as.POSIXct(data_ini - 1, tz = "UTC") &
               dwd_h$hora <  as.POSIXct(data_fim + 2, tz = "UTC"), ]
dwd_d <- diario(dwd_h)

# ---- bruce data contract ----------------------------------------------------------
diarios <- rbind(cbind(fair_d, source = "fair"), cbind(dwd_d, source = "dwd"))
diarios <- diarios[diarios$date >= data_ini & diarios$date <= data_fim, ]
obs <- data.frame(station_id = diarios$station_id, date = diarios$date, value = diarios[[VAR]])
stations <- rbind(
  data.frame(station_id = fair_st$station_id, source = "fair", lon = fair_st$lon, lat = fair_st$lat),
  data.frame(station_id = dwd_sel$station_id, source = "dwd", lon = dwd_sel$lon, lat = dwd_sel$lat))
stations <- stations[stations$station_id %in% obs$station_id, ]
sources <- data.frame(source = c("dwd", "fair"), bias = c(FALSE, TRUE), obs_sd = c(0.3, NA))

# ---- ERA5-Land hourly -> daily on local days ----------------------------------------
bb <- st_bbox(st_transform(area, 4326))
nc <- file.path(DIR_OUT, "cache", sprintf("era5land_t2m_%d-%02d.nc", ANO, MES))
if (!file.exists(nc)) {
  dias <- format(seq(data_ini - 1, data_fim + 1, by = "day"))   # one day margin for local days
  req <- list(dataset_short_name = "reanalysis-era5-land", variable = "2m_temperature",
              year = unique(substr(dias, 1, 4)), month = unique(substr(dias, 6, 7)),
              day = sprintf("%02d", 1:31), time = sprintf("%02d:00", 0:23),
              area = round(c(bb["ymax"], bb["xmin"], bb["ymin"], bb["xmax"]) + c(.1, -.1, -.1, .1), 2),
              data_format = "netcdf", download_format = "unarchived", target = basename(nc))
  ecmwfr::wf_request(req, path = dirname(nc))
}
e5 <- rast(nc) - 273.15
t_local <- as.Date(format(time(e5), tz = TZ, usetz = FALSE))
fun <- switch(VAR, median = stats::median, max = max, min = min)
e5d <- tapp(e5, index = as.character(t_local), fun = fun)
nm <- gsub("[^0-9]", "-", sub("^X", "", names(e5d)))
dias_e5 <- as.Date(nm, format = "%Y-%m-%d")
keep <- dias_e5 >= data_ini & dias_e5 <= data_fim
xy <- as.data.frame(e5d[[which(keep)]], xy = TRUE, na.rm = TRUE)
background <- do.call(rbind, lapply(seq_len(sum(keep)), function(k)
  data.frame(lon = xy$x, lat = xy$y, date = dias_e5[keep][k], value = xy[[2 + k]])))
inside <- st_intersects(st_transform(st_as_sf(background, coords = c("lon", "lat"), crs = 4326), CRS_UTM),
                        area, sparse = FALSE)[, 1]
background <- background[inside, ]

# ==================================================================
# 2. Covariates available anywhere
# ==================================================================
dem <- elevatr::get_elev_raster(locations = st_transform(area, 4326), z = 10, clip = "locations") |>
  rast() |> project(paste0("EPSG:", CRS_UTM))

wc <- rast("/vsicurl/https://esa-worldcover.s3.eu-central-1.amazonaws.com/v200/2021/map/ESA_WorldCover_10m_2021_v200_N51E012_Map.tif")
wc <- crop(wc, ext(bb["xmin"], bb["xmax"], bb["ymin"], bb["ymax"]))
urb <- aggregate(wc == 50, fact = 50, fun = "mean", na.rm = TRUE)   # built-up fraction, ~500 m
urb <- project(urb, paste0("EPSG:", CRS_UTM), method = "average")
names(urb) <- "urb"

lcz <- LCZ4r::lcz_get_map(city = "Berlin")
svf <- LCZ4r::lcz_get_parameters(lcz, istack = FALSE, iselect = "SVFmean")

# ==================================================================
# 3. bruce: same steps as for Rio, with the generic defaults
# ==================================================================
d <- bruce_data(obs, stations, sources, crs = CRS_UTM, background = background, background_sd = 0.7)
d <- bruce_covariates(d, list(elev = dem, urb = urb, svf = svf),
                      background_method = c(elev = "buffer_mean", urb = "buffer_mean", svf = "zero"),
                      background_buffer = 5500, na_value = c(urb = 0))
d

mesh <- bruce_mesh(d)                       # sizes derived from the station spacing
spec <- bruce_spec()                        # range prior anchored on the network ("nn")
fit  <- bruce_fit(d, mesh, spec, keep_config = TRUE)
summary(fit)

grid500 <- rast(vect(berlim), resolution = 500, crs = paste0("EPSG:", CRS_UTM))
values(grid500) <- 1
grid500 <- mask(grid500, vect(berlim))
rs <- c(elev = "average", urb = "near", svf = "near")
pred <- bruce_predict(fit, grid = grid500, n_samples = 200, resample = rs)
dec  <- bruce_decompose(fit, grid = grid500, resample = rs)
print(dec)

dia_quente <- names(which.max(tapply(obs$value, obs$date, mean)))
png(file.path(DIR_OUT, sprintf("decomposicao_%s.png", tag)), width = 1500, height = 520, res = 110)
plot(dec, date = dia_quente)
dev.off()

# Leave-station-out on the DWD reference stations
cv <- bruce_cv(d, mesh, spec, scheme = "lso", n_samples = 200,
               checkpoint_dir = file.path(DIR_OUT, sprintf("cv_lso_%s", tag)))
print(cv)

saveRDS(list(fit = fit, pred_mean = wrap(pred$mean), pred_sd = wrap(pred$sd),
             variance = dec$variance, cv = cv),
        file.path(DIR_OUT, sprintf("bruce_berlin_%s.rds", tag)))
writeRaster(pred$mean, file.path(DIR_OUT, sprintf("bruce_berlin_pred_%s.tif", tag)), overwrite = TRUE)
message("Done: ", DIR_OUT)
