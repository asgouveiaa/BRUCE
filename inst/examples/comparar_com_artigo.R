# ==================================================================
# comparar_com_artigo.R
# Compares the bruce fit from rio_adapter.R with the article's saved fit
# for the same month (modelo_9meses_median_max.R output). Run from
# ~/scratch/spde_paper after rio_adapter.R.
#
# What "the same" means here: INLA results are not bit-identical between
# runs and the adapter differs in a few small details (see its header),
# so we expect
#   * fixed effects within ~0.5 posterior sd of the original,
#   * hyperparameters of the same order (Range, Stdev, GroupRho),
#   * daily maps with correlation > 0.99 and RMSE of a few tenths of a degree.
# ==================================================================

library(terra)

VAR <- "median"; MES <- "2025-01"
tag <- sprintf("%s_%s", VAR, MES)

orig  <- readRDS(sprintf("resultados_9meses_median_max/pacote_mensal_independente_%s_%s.rds", VAR, MES))
novo  <- readRDS(sprintf("resultados_bruce/bruce_%s.rds", tag))
fit_o <- orig$fit
fit_n <- novo$fit$inla

# ---- 1. Fixed effects ----------------------------------------------------
nomes <- c(b0 = "intercept", elev_std = "elev", urb_std = "urb", svf_std = "svf",
           aura = "bias_aurassure", wunder = "bias_wunder", aeroporto = "bias_airport")
nomes <- nomes[names(nomes) %in% rownames(fit_o$summary.fixed) &
               nomes %in% rownames(fit_n$summary.fixed)]
fx <- data.frame(
  termo       = names(nomes),
  artigo      = fit_o$summary.fixed[names(nomes), "mean"],
  bruce       = fit_n$summary.fixed[nomes, "mean"],
  sd_artigo   = fit_o$summary.fixed[names(nomes), "sd"]
)
fx$diff_em_sd <- (fx$bruce - fx$artigo) / fx$sd_artigo
cat("\n=== Efeitos fixos ===\n"); print(fx, digits = 3, row.names = FALSE)

# ---- 2. Hyperparameters --------------------------------------------------
hp <- function(f, padrao) {
  i <- grep(padrao, rownames(f$summary.hyperpar))
  if (length(i)) f$summary.hyperpar[i[1], "mean"] else NA
}
hps <- data.frame(
  parametro = c("Range (m)", "Stdev campo", "GroupRho", "Rho diario", "Prec diario", "Prec estacao"),
  artigo = c(hp(fit_o, "^Range for s"), hp(fit_o, "^Stdev for s"), hp(fit_o, "^GroupRho for s"),
             hp(fit_o, "^Rho for dia_idx"), hp(fit_o, "^Precision for dia_idx"),
             hp(fit_o, "^Precision for id_estacao_idx")),
  bruce  = c(hp(fit_n, "^Range for s"), hp(fit_n, "^Stdev for s"), hp(fit_n, "^GroupRho for s"),
             hp(fit_n, "^Rho for day_idx"), hp(fit_n, "^Precision for day_idx"),
             hp(fit_n, "^Precision for station_idx"))
)
cat("\n=== Hiperparametros ===\n"); print(hps, digits = 4, row.names = FALSE)
cat("\nStatus do ajuste bruce:", novo$fit$status, "\n")

# ---- 3. Daily maps -------------------------------------------------------
pred_o <- read.csv(sprintf("resultados_9meses_median_max/predicao_mensal_independente_%s_%s.csv", VAR, MES))
pred_o$data <- as.Date(pred_o$data)
pred_n <- unwrap(novo$pred_mean)

mapas <- do.call(rbind, lapply(sort(unique(pred_o$data)), function(dia) {
  o <- pred_o[pred_o$data == dia, ]
  n <- terra::extract(pred_n[[as.character(dia)]], cbind(o$x, o$y))[, 1]
  ok <- !is.na(n)
  data.frame(data = dia, n_pixels = sum(ok),
             cor = cor(o$pred[ok], n[ok]),
             rmse = sqrt(mean((n[ok] - o$pred[ok])^2)),
             vies = mean(n[ok] - o$pred[ok]),
             max_abs = max(abs(n[ok] - o$pred[ok])))
}))
cat("\n=== Mapas diarios (bruce - artigo) ===\n"); print(mapas, digits = 3, row.names = FALSE)
cat(sprintf("\nResumo: correlacao media %.4f | RMSE medio %.3f C | maior diferenca %.2f C\n",
            mean(mapas$cor), mean(mapas$rmse), max(mapas$max_abs)))

# ---- 4. Map of the difference on the worst day ----------------------------
pior <- mapas$data[which.max(mapas$rmse)]
o <- pred_o[pred_o$data == pior, ]
r_o <- rast(o[, c("x", "y", "pred")], type = "xyz", crs = crs(pred_n))
r_n <- resample(pred_n[[as.character(pior)]], r_o, method = "near")
png(sprintf("resultados_bruce/comparacao_%s.png", tag), width = 1500, height = 520, res = 110)
op <- par(mfrow = c(1, 3), mar = c(2, 2, 3, 5))
rng <- range(c(values(r_o), values(r_n)), na.rm = TRUE)
plot(r_o, main = sprintf("Artigo (%s)", pior), range = rng, col = hcl.colors(50, "YlOrRd", rev = TRUE))
plot(r_n, main = "bruce", range = rng, col = hcl.colors(50, "YlOrRd", rev = TRUE))
dd <- r_n - r_o; lim <- max(abs(values(dd)), na.rm = TRUE)
plot(dd, main = "bruce - artigo", range = c(-lim, lim), col = hcl.colors(50, "Blue-Red 3"))
par(op); dev.off()
cat("Figura:", sprintf("resultados_bruce/comparacao_%s.png", tag), "\n")
