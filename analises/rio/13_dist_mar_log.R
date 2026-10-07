# ==================================================================
# 13_dist_mar_log.R
# Distance to the sea: straight line or log(1 + dist)? Fits nothing; takes
# seconds. Uses the distance already in the cache (sea = DEM <= 0 m).
#
# The log makes the first kilometres count more than the last ones: going
# from 0 to 1 km from the sea changes much more than going from 20 to 21 km.
#
# Quick check without INLA: the mean Td of each station against
# network + elevation + distance, in the two versions. Lower residual sd and
# AIC = the distance explains the stations better. The final word is the
# full model: 12_orvalho_ajuste.R dist_log=1.
#
# Needs 12_orvalho_ajuste.R so_preparar=1 already run (cache).
# Rode de dentro de analises/rio:  Rscript 13_dist_mar_log.R
# ==================================================================

library(terra)
source("config.R")

CRS_UTM <- 31983
DIR_CACHE <- file.path(RES, "orvalho", "cache")
VAR <- "median"; ANO <- 2026; MES <- 1

dist_mar <- rast(file.path(DIR_CACHE, "dist_mar.tif"))
dem <- rast(file.path(DIR_CACHE, "dem.tif"))

td <- readRDS(file.path(DIR_CACHE, sprintf("td_obs_%s_%d-%02d.rds", VAR, ANO, MES)))
col <- paste0(VAR, "_orvalho")
td <- td[!is.na(td$latitude) & !is.na(td$longitude) & !is.na(td[[col]]), ]
td$id <- paste(td$origem, td$estacao_id, sep = ":")
td$val <- td[[col]]
est <- aggregate(val ~ id + origem + longitude + latitude, data = td, FUN = mean)
names(est)[names(est) == "val"] <- "td"
est <- est[!duplicated(est$id), ]
pts <- project(vect(est, geom = c("longitude", "latitude"), crs = "EPSG:4326"), paste0("EPSG:", CRS_UTM))
est$dist <- extract(dist_mar, pts)[, 2]
est$log_dist <- log1p(est$dist)
est$elev <- extract(dem, pts)[, 2]

modelos <- list(reta = lm(td ~ origem + elev + dist, data = est),
                log  = lm(td ~ origem + elev + log_dist, data = est))
tab <- do.call(rbind, lapply(names(modelos), function(nm) {
  m <- modelos[[nm]]; co <- summary(m)$coefficients
  v <- c(reta = "dist", log = "log_dist")[[nm]]
  if (!v %in% rownames(co)) stop("A distancia ficou colinear com as outras variaveis")
  data.frame(versao = nm, coef_dist = round(co[v, 1], 3), t = round(co[v, 3], 1),
             sd_residuo = round(sigma(m), 3), r2_aj = round(summary(m)$adj.r.squared, 3),
             AIC = round(AIC(m), 1))
}))
cat("=== Td medio da estacao ~ rede + elevacao + distancia ao mar ===\n")
print(tab, row.names = FALSE)
cat("Menor sd_residuo / AIC = explica melhor. Diferenca de AIC < 2: empate.\n")
cat("(So", nrow(est), "estacoes e sem campo espacial: checagem rapida, nao o modelo.)\n")

cat("\n=== INMET e aeroportos: residuo (observado - previsto, C) ===\n")
ref <- est$origem %in% c("inmet", "redemet")
print(data.frame(id = est$id[ref], dist_km = round(est$dist[ref], 2), td = round(est$td[ref], 1),
                 res_reta = round(residuals(modelos$reta)[ref], 2),
                 res_log = round(residuals(modelos$log)[ref], 2)), row.names = FALSE)

png("dist_mar_reta_log.png", width = 1200, height = 520, res = 110)
par(mfrow = c(1, 2))
parcial <- function(m, x, nome) {
  r <- residuals(m) + coef(m)[nome] * est[[nome]]
  plot(est[[x]], r, pch = 20, col = ifelse(ref, "red", "grey40"), xlab = "distancia ao mar (km)",
       ylab = "Td parcial (C)", main = sprintf("%s (vermelho = INMET/aeroporto)", nome))
  o <- order(est[[x]]); lines(est[[x]][o], (coef(m)[nome] * est[[nome]])[o], lwd = 2)
}
parcial(modelos$reta, "dist", "dist")
parcial(modelos$log, "dist", "log_dist")
dev.off()
cat("\nFigura: dist_mar_reta_log.png (a curva mostra a forma do efeito em cada versao)\n")
