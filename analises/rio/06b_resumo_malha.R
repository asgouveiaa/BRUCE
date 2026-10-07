# ==================================================================
# 06b_resumo_malha.R
# Summary table and figure of 06_sensibilidade_malha.R from the meshes that
# have finished so far (reads the saved malha_*.rds; refits nothing).
# Can run while 06 is still running. Same settings as the top of 06.
# Rode de dentro de analises/rio:  Rscript 06b_resumo_malha.R
# ==================================================================

library(terra)
source("config.R")

VAR <- "median"; MES <- "2025-01"
INT_STRATEGY <- "eb"
N_DIAS       <- NULL
DIR <- file.path(RES, sprintf("sensibilidade_malha_%s%s", INT_STRATEGY,
                              if (is.null(N_DIAS)) "" else sprintf("_%dd", N_DIAS)))
if (!length(list.files(DIR, "^malha_.*\\.rds$"))) stop("Nenhuma malha terminada em ", DIR)

pred_art <- read.csv(sprintf(file.path(ARTIGO, "predicao_mensal_independente_%s_%s.csv"), VAR, MES))
pred_art$data <- as.Date(pred_art$data)

# ==================================================================
# Summary
# ==================================================================
res <- lapply(list.files(DIR, "^malha_[^_]*\\.rds$", full.names = TRUE), readRDS)   # skips *_colapsada.rds
res <- res[order(vapply(res, `[[`, 0, "n"))]
fina <- unwrap(res[[length(res)]]$pred_mean)

comparar <- function(r, ref) {
  dias <- sort(unique(pred_art$data))
  dias <- dias[as.character(dias) %in% names(unwrap(r$pred_mean))]
  cmp <- vapply(dias, function(dia) {
    p <- pred_art[pred_art$data == dia, ]
    a <- terra::extract(unwrap(r$pred_mean)[[as.character(dia)]], cbind(p$x, p$y))[, 1]
    b <- terra::extract(ref[[as.character(dia)]], cbind(p$x, p$y))[, 1]
    c(cor_art = cor(a, p$pred, use = "complete"), rmse_art = sqrt(mean((a - p$pred)^2, na.rm = TRUE)),
      rmse_fina = sqrt(mean((a - b)^2, na.rm = TRUE)))
  }, numeric(3))
  rowMeans(cmp)
}
tab <- do.call(rbind, lapply(res, function(r) {
  h <- r$hyper; m <- comparar(r, fina)
  data.frame(malha = r$nome, status = r$status, vertices = r$n, minutos = round(r$minutos, 1),
             range = h["Range for s", "mean"], range_q025 = h["Range for s", "0.025quant"],
             range_q975 = h["Range for s", "0.975quant"],
             stdev_campo = h["Stdev for s", "mean"], grouprho = h["GroupRho for s", "mean"],
             aresta_tipica = r$mesh_resolution$edge, range_sobre_aresta = r$mesh_resolution$ratio,
             elev = r$fixed["elev", "mean"], cor_mapa_artigo = m["cor_art"],
             rmse_mapa_artigo = m["rmse_art"], rmse_mapa_vs_mais_fina = m["rmse_fina"])
}))
rownames(tab) <- NULL
cat("\n=== Sensibilidade a malha ===\n"); print(tab, digits = 3)
write.csv(tab, file.path(DIR, "resumo_sensibilidade_malha.csv"), row.names = FALSE)

png(file.path(DIR, "range_vs_malha.png"), width = 900, height = 600, res = 110)
par(mar = c(4.5, 4.5, 1, 1))
plot(tab$vertices, tab$range / 1000, log = "x", pch = 19, ylim = range(c(tab$range_q025, tab$range_q975)) / 1000,
     xlab = "Vertices da malha (escala log)", ylab = "Range posterior (km), media e IC 95%")
arrows(tab$vertices, tab$range_q025 / 1000, tab$vertices, tab$range_q975 / 1000, angle = 90, code = 3, length = 0.04)
text(tab$vertices, tab$range_q975 / 1000, tab$malha, pos = 3, cex = 0.8)
dev.off()
cat("Figura:", file.path(DIR, "range_vs_malha.png"), "\n")
