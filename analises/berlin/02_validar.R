# ==================================================================
# 02_validar.R
# Does the bruce map of Berlin add anything over simply using the
# official DWD stations?
#
# Uses the leave-station-out already run by 01_ajuste_mes.R (each DWD
# station left out in turn, the model refitted without it, and its
# temperature predicted). The same held-out station-days are predicted by
# two simple alternatives that use only the other DWD stations:
#   media_dwd     mean of the other DWD stations on that day
#   vizinha_dwd   the nearest other DWD station on that day
# If bruce does not beat these, the map is not adding information here.
#
# Also checks the uncertainty: the 90% intervals of the cross-validation
# describe the spatial field only; adding the station's own offset and the
# measurement error gives the interval for what a station would read.
#
# Rode de dentro de analises/berlin (resultados de 01 em resultados/):
#   cd ~/bruce/analises/berlin && Rscript 02_validar.R
# Takes seconds (no refitting). Results in resultados/validacao_<tag>/.
# ==================================================================

library(terra)

ANO <- 2022; MES <- 7; VAR <- "median"
tag <- sprintf("%s_%d-%02d", VAR, ANO, MES)
DIR_RES <- "resultados"
DIR_OUT <- file.path(DIR_RES, sprintf("validacao_%s", tag))
dir.create(DIR_OUT, showWarnings = FALSE, recursive = TRUE)

arq <- file.path(DIR_RES, sprintf("bruce_berlin_%s.rds", tag))
if (!file.exists(arq)) stop("Nao achei ", arq, ". Rode de dentro de analises/berlin.")
cat("Lendo", arq, "(pode levar alguns segundos)\n")
r  <- readRDS(arq)
d  <- r$fit$data
cv <- r$cv
pr <- as.data.frame(cv$predictions)

# ---- 1. Simple alternatives on the same held-out station-days ---------------
dwd <- d$obs[d$obs$source == "dwd", c("station_id", "date", "value", "x", "y")]
alt <- t(vapply(seq_len(nrow(pr)), function(i) {
  o <- dwd[dwd$date == pr$date[i] & dwd$station_id != pr$station_id[i], ]
  if (!nrow(o)) return(c(NA_real_, NA_real_))
  me <- dwd[dwd$station_id == pr$station_id[i], c("x", "y")][1, ]
  dist <- sqrt((o$x - me$x)^2 + (o$y - me$y)^2)
  c(mean(o$value), o$value[which.min(dist)])
}, numeric(2)))
pr$media_dwd   <- alt[, 1]
pr$vizinha_dwd <- alt[, 2]
pr <- pr[stats::complete.cases(pr[, c("pred", "media_dwd", "vizinha_dwd")]), ]

erro <- function(p, o) c(rmse = sqrt(mean((p - o)^2)), mae = mean(abs(p - o)), vies = mean(p - o))
tab <- rbind(bruce       = erro(pr$pred, pr$obs),
             media_dwd   = erro(pr$media_dwd, pr$obs),
             vizinha_dwd = erro(pr$vizinha_dwd, pr$obs))
tab <- data.frame(metodo = rownames(tab), round(tab, 3), row.names = NULL)
tab$ganho_bruce_pct <- round(100 * (1 - tab$rmse[1] / tab$rmse), 1)   # how much smaller bruce's RMSE is
cat("\n=== 1. bruce contra alternativas simples (", nrow(pr), " estacao-dias) ===\n", sep = "")
print(tab, row.names = FALSE)
cat("ganho_bruce_pct: quanto o RMSE do bruce e menor que o do metodo (positivo = bruce melhor)\n")
write.csv(tab, file.path(DIR_OUT, "comparacao_metodos.csv"), row.names = FALSE)

# ---- 2. Per station: is the error spread or concentrated? -------------------
por_est <- do.call(rbind, lapply(split(pr, pr$station_id), function(p) data.frame(
  estacao = p$station_id[1], n = nrow(p),
  rmse_bruce = sqrt(mean((p$pred - p$obs)^2)), vies_bruce = mean(p$pred - p$obs),
  rmse_media_dwd = sqrt(mean((p$media_dwd - p$obs)^2)),
  rmse_vizinha_dwd = sqrt(mean((p$vizinha_dwd - p$obs)^2)))))
por_est[, -(1:2)] <- round(por_est[, -(1:2)], 3)
por_est$bruce_melhor <- por_est$rmse_bruce < pmin(por_est$rmse_media_dwd, por_est$rmse_vizinha_dwd)
por_est <- por_est[order(-abs(por_est$vies_bruce)), ]
cat("\n=== 2. Por estacao (ordenado pelo vies do bruce) ===\n")
print(por_est, row.names = FALSE)
cat(sprintf("bruce melhor que as duas alternativas em %d de %d estacoes\n",
            sum(por_est$bruce_melhor), nrow(por_est)))
write.csv(por_est, file.path(DIR_OUT, "por_estacao.csv"), row.names = FALSE)

# ---- 3. Coverage of the 90% interval ------------------------------------------
hp <- r$fit$inla$summary.hyperpar
sd_estacao <- 1 / sqrt(hp["Precision for station_idx", "mean"])
sd_medida  <- d$sources$obs_sd[d$sources$source == "dwd"]
z <- qnorm(1 - (1 - cv$level) / 2)
sd_total <- sqrt(pr$pred_sd^2 + sd_estacao^2 + sd_medida^2)
cob <- data.frame(
  intervalo = c("so o campo (como saiu do bruce_cv)", "campo + jeito da estacao + erro de medida"),
  cobertura = round(c(mean(pr$obs >= pr$lower & pr$obs <= pr$upper),
                      mean(abs(pr$obs - pr$pred) <= z * sd_total)), 3),
  largura_media_C = round(c(mean(pr$upper - pr$lower), mean(2 * z * sd_total)), 2))
cat(sprintf("\n=== 3. Cobertura do intervalo de %.0f%% (o ideal e perto de %.2f) ===\n",
            100 * cv$level, cv$level))
print(cob, row.names = FALSE)
cat(sprintf("(sd do jeito da estacao %.2f C, erro de medida da DWD %.2f C)\n", sd_estacao, sd_medida))
write.csv(cob, file.path(DIR_OUT, "cobertura.csv"), row.names = FALSE)

# ---- 4. Figures ----------------------------------------------------------------
# Observed vs predicted, the three methods side by side
png(file.path(DIR_OUT, "observado_vs_predito.png"), width = 1500, height = 520, res = 110)
par(mfrow = c(1, 3), mar = c(4.5, 4.5, 2.5, 1))
lim <- range(pr[, c("obs", "pred", "media_dwd", "vizinha_dwd")])
for (m in c("pred", "media_dwd", "vizinha_dwd")) {
  plot(pr$obs, pr[[m]], xlim = lim, ylim = lim, pch = 19, cex = 0.5, col = "#00000060",
       xlab = "Observado na estacao deixada de fora (C)", ylab = "Predito (C)",
       main = sprintf("%s  (RMSE %.2f)", c(pred = "bruce", media_dwd = "media das outras DWD",
                                             vizinha_dwd = "DWD mais proxima")[m],
                      sqrt(mean((pr[[m]] - pr$obs)^2))))
  abline(0, 1, col = "red")
}
dev.off()

# Every day as a difference from that day's city mean, all on the same scale
p <- rast(file.path(DIR_RES, sprintf("bruce_berlin_pred_%s.tif", tag)))
anom <- p - global(p, "mean", na.rm = TRUE)$mean
lim_a <- max(abs(global(anom, "range", na.rm = TRUE)))
nd <- nlyr(anom)
png(file.path(DIR_OUT, "anomalia_diaria.png"), width = 2100, height = 330 * ceiling(nd / 7), res = 110)
par(mfrow = c(ceiling(nd / 7), 7), mar = c(0.5, 0.5, 1.8, 3.5))
for (k in seq_len(nd)) {
  plot(anom[[k]], range = c(-lim_a, lim_a), col = hcl.colors(50, "Blue-Red 3"), axes = FALSE,
       main = names(anom)[k], legend = (k %% 7 == 0 || k == nd))   # one colour bar per row
}
dev.off()
amp <- global(p, "range", na.rm = TRUE)
cat(sprintf("\nDiferenca entre o ponto mais quente e o mais frio da cidade: de %.1f a %.1f C (media %.1f C)\n",
            min(amp[, 2] - amp[, 1]), max(amp[, 2] - amp[, 1]), mean(amp[, 2] - amp[, 1])))
cat(sprintf("Escala comum dos mapas de anomalia: -%.1f a +%.1f C em relacao a media da cidade no dia\n",
            lim_a, lim_a))
cat("\nArquivos em", DIR_OUT, "\n")
