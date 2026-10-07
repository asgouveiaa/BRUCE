# ==================================================================
# 15_orvalho_efeito_estacao.R
# Is the dew point map too smooth? Is a covariate missing?
#
# The model gives each station its own fixed offset ("station effect": the
# station reads a bit more or less humid than the map, every day of the
# month). That offset is a mix of two things the model cannot tell apart:
#   - the sensor (calibration, shelter, installation), and
#   - REAL local climate the covariates do not explain (next to a forest,
#     a lagoon, a big paved area).
# If part of it is real climate, the offsets should line up with something on
# the ground. Here they are compared with candidate covariates measured
# around each station (500 m and 1 km): vegetation, water, urban fraction,
# plus the ones already in the model (as a check: they should show nothing).
#
# Two more tests need no covariate at all: whether the official stations
# (good sensors) have offsets as large as the citizen ones, and whether
# NEIGHBOURING stations of different networks have similar offsets (different
# sensors cannot share a calibration error, so agreement = real climate).
#
# Uses the fit already saved by 12 (nothing is refitted; about a minute,
# most of it reading the land use shapefile).
# Rode de dentro de analises/rio:
#   Rscript 15_orvalho_efeito_estacao.R
#   Rscript 15_orvalho_efeito_estacao.R tag=td_median_2023-11_buf10km_grho0.9
# ==================================================================

library(sf)
library(terra)
library(bruce)
source("config.R")

args <- commandArgs(trailingOnly = TRUE)
arg <- function(nome, padrao) {
  a <- grep(paste0("^", nome, "="), args, value = TRUE)
  if (length(a)) sub("^[^=]*=", "", a[1]) else padrao
}
TAG <- arg("tag", "td_median_2026-01_buf10km")
CRS_UTM <- 31983
DIR_OUT <- file.path(RES, "orvalho")
RAIOS <- c(500, 1000)
# land use classes (column usoagregad) taken as vegetation / water; the
# classes found are printed, so this can be corrected
PADRAO_VEG  <- "arb\u00f3re|arbore|arbust|florest|veget|gram\u00edn|gramin|mangue|restinga|campo|agr\u00edc|agric"
PADRAO_AGUA <- "água|agua|hídric|hidric|lagoa|corpo"
# urban = the same classes used in 12 (matched without accents, so encoding
# differences in the shapefile do not matter)
PADRAO_URB <- "residenc|comerc|industr|transporte|institucion|educa|favela"
sem_acento <- function(x) iconv(enc2utf8(x), "UTF-8", "ASCII//TRANSLIT", sub = "")

fit <- readRDS(file.path(DIR_OUT, sprintf("bruce_%s.rds", TAG)))$fit

# ---- station effects --------------------------------------------------------------
re <- fit$inla$summary.random$station_idx
ef <- data.frame(station_id = fit$station_ids[re$ID], efeito = re$mean, efeito_sd = re$sd)
est <- merge(fit$data$stations, ef, by = "station_id")
sd_ef <- 1 / sqrt(fit$inla$summary.hyperpar["Precision for station_idx", "mean"])
cat(sprintf("Ajuste %s: %d estacoes; sd do efeito de estacao = %.2f C\n", TAG, nrow(est), sd_ef))
cat("Efeito medio por rede (o vies da rede ja foi tirado; deve dar ~0):\n")
print(round(tapply(est$efeito, est$source, mean), 2))

# ---- candidate covariates around each station -------------------------------------------
pts <- st_as_sf(est, coords = c("x", "y"), crs = CRS_UTM, remove = FALSE)
uso <- st_transform(read_sf(USO_SOLO), CRS_UTM)
classes <- sort(unique(uso$usoagregad))
rot <- paste(uso$usoagregad, sem_acento(uso$usoagregad))
uso$veg  <- grepl(PADRAO_VEG, rot, ignore.case = TRUE)
uso$agua <- grepl(PADRAO_AGUA, rot, ignore.case = TRUE)
uso$urb  <- grepl(PADRAO_URB, rot, ignore.case = TRUE)
cat("\nClasses do uso do solo (V = vegetacao, A = agua, U = urbano):\n")
for (cl in classes) {
  u <- uso[uso$usoagregad == cl, ][1, ]
  cat(sprintf("  [%s%s%s] %s\n", if (u$veg) "V" else " ", if (u$agua) "A" else " ",
              if (u$urb) "U" else " ", cl))
}
uso <- uso[uso$veg | uso$agua | uso$urb, c("veg", "agua", "urb")]
uso <- st_make_valid(uso)

fracao <- function(raio) {
  buf <- st_buffer(pts, raio)
  inter <- suppressWarnings(st_intersection(buf[, "station_id"], uso))
  inter$area <- as.numeric(st_area(inter))
  a_tot <- pi * raio^2
  agg <- function(col) tapply(inter$area * inter[[col]], inter$station_id, sum) / a_tot
  out <- data.frame(station_id = pts$station_id)
  for (col in c("veg", "agua", "urb")) {
    v <- agg(col); out[[paste0(col, "_", raio)]] <- round(as.numeric(v[out$station_id]), 3)
  }
  out[is.na(out)] <- 0
  out
}
for (r in RAIOS) est <- merge(est, fracao(r), by = "station_id")

arq_cov <- file.path(DIR_OUT, "cache", c(elev = "dem.tif", dist_mar = "dist_mar.tif"))
v <- vect(est, geom = c("x", "y"), crs = paste0("EPSG:", CRS_UTM))
est$elev <- extract(rast(arq_cov[1]), v)[, 2]
est$dist_mar <- extract(rast(arq_cov[2]), v)[, 2]

# ---- does anything line up with the station effects? ------------------------------------
cands <- c(paste0(c("veg", "agua", "urb"), "_", rep(RAIOS, each = 3)), "elev", "dist_mar")
tab <- do.call(rbind, lapply(cands, function(cn) {
  if (sd(est[[cn]], na.rm = TRUE) == 0) return(NULL)
  m <- lm(as.formula(paste("efeito ~ source +", cn)), data = est, weights = 1 / est$efeito_sd^2)
  co <- summary(m)$coefficients[cn, ]
  data.frame(candidata = cn, faixa = sprintf("%.2f-%.2f", min(est[[cn]], na.rm = TRUE), max(est[[cn]], na.rm = TRUE)),
             cor = round(cor(est$efeito, est[[cn]], use = "complete.obs"), 2),
             coef = round(co[1], 3), t = round(co[3], 1), p = signif(co[4], 2))
}))
cat("\n=== Efeito da estacao ~ rede + candidata (cada uma sozinha) ===\n")
print(tab, row.names = FALSE)
cat("elev e dist_mar ja estao no modelo: devem dar ~0 (checagem).\n")
cat("Uma candidata nova com |t| > 2 sugere clima local que o mapa nao capta.\n")
cat(sprintf("Correlacao vegetacao (1 km) x elevacao nas estacoes: %.2f\n",
            cor(est$veg_1000, est$elev, use = "complete.obs")))

cat("\n=== Estacoes com maior efeito (lidas mais umidas (+) ou mais secas (-) que o mapa) ===\n")
print(head(est[order(-abs(est$efeito)), c("station_id", "source", "efeito", "veg_1000", "agua_1000",
                                          "urb_1000", "elev")], 12), row.names = FALSE)
# ---- sensor or real climate? two tests that need no covariate -------------------------
# (a) Official stations (INMET, airports) have good, calibrated sensors: their
#     offsets are mostly REAL local climate. If they are as large as the
#     citizen stations' offsets, real climate is a big part of the offsets.
cat("\n=== (a) Tamanho do efeito por rede (desvio-padrao dos efeitos, C) ===\n")
print(round(tapply(est$efeito, est$source, function(v) sqrt(mean(v^2))), 2))
cat("INMET/aeroportos tem sensor bom: o efeito deles e quase so clima local real.\n")

# (b) Sensor defects are independent from one station to the next. Real local
#     climate is shared by neighbours. So: do nearby stations have similar
#     offsets, above all nearby stations of DIFFERENT networks (different
#     sensors cannot share a calibration error)?
D <- as.matrix(dist(est[, c("x", "y")]))
viz <- function(ate_m, so_redes_diferentes) {
  w <- (D > 0 & D <= ate_m) * 1
  if (so_redes_diferentes) w <- w * outer(est$source, est$source, "!=")
  w
}
moran <- function(z, w) {
  z <- z - mean(z); sum(w * outer(z, z)) / sum(w) / mean(z^2)
}
set.seed(1)
cat("\n=== (b) Vizinhos tem efeitos parecidos? (I de Moran; acaso = ~0) ===\n")
for (r in c(2000, 4000)) for (dif in c(FALSE, TRUE)) {
  w <- viz(r, dif)
  if (sum(w) == 0) next
  obs_i <- moran(est$efeito, w)
  perm <- replicate(999, moran(sample(est$efeito), w))
  cat(sprintf("  ate %d km, %-26s pares %3d   I = %5.2f   p = %.3f\n", r / 1000,
              if (dif) "so redes diferentes:" else "todas as redes:", sum(w) / 2, obs_i,
              mean(c(perm, obs_i) >= obs_i)))
}
cat("I positivo com p < 0.05 em 'so redes diferentes' = sensores diferentes concordam:\n")
cat("isso e clima local real que o mapa esta deixando no efeito da estacao.\n")

write.csv(est, file.path(DIR_OUT, sprintf("efeito_estacao_%s.csv", TAG)), row.names = FALSE)

# ---- map ---------------------------------------------------------------------------
png(file.path(DIR_OUT, sprintf("efeito_estacao_%s.png", TAG)), width = 1300, height = 650, res = 110)
pal <- colorRampPalette(c("#b2182b", "#f7f7f7", "#2166ac"))(21)
lim <- max(abs(est$efeito))
cor_pt <- pal[cut(est$efeito, seq(-lim, lim, length.out = 22), include.lowest = TRUE)]
titulo <- "Efeito de cada estacao (azul = mais umida que o mapa, vermelho = mais seca)"
mrj <- tryCatch(st_transform(geobr::read_municipality(3304557, year = 2020, showProgress = FALSE), CRS_UTM),
                error = function(e) NULL)
if (is.null(mrj)) plot(est$x, est$y, type = "n", asp = 1, xlab = "", ylab = "", main = titulo) else
  plot(st_geometry(mrj), border = "grey50", main = titulo)
points(est$x, est$y, pch = c(airport = 21, inmet = 23, alerta_rio = 24, alertario = 24, aurassure = 22,
                             wunder = 25)[est$source], bg = cor_pt, cex = 1.6)
legend("bottomleft", c("aeroporto", "INMET", "Alerta Rio", "Aurassure", "Wunderground"),
       pch = c(21, 23, 24, 22, 25), bty = "n", cex = 0.8)
dev.off()
cat("\nSalvos:", sprintf("efeito_estacao_%s.csv e .png", TAG), "em", DIR_OUT, "\n")
