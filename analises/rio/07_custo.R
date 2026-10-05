# ==================================================================
# 07_custo.R
# Where does the cost of bruce go, and what can be cut without
# changing the model? January 2025 median, same data, mesh and priors
# as 01_ajuste_mes.R (read from the fit it saved).
#
# Parte A - one fit, several ways to run INLA:
#   referencia      as in 01 (ccd, INLA's default threads)
#   threads_*       other thread layouts (num.threads = "A:B")
#   pardiso         PARDISO solver, if licensed on this machine
#   eb              int.strategy = "eb" (posterior mode only)
#   eb_gaussian     eb + Gaussian approximation of the latent marginals
#   theta_fixo      hyperparameters fixed at the mode of the 01 fit (the
#                   cost of a daily operational map with theta learned
#                   once a month)
#   Each is compared with the reference: time, hyperparameters, fixed
#   effects, field sd, and the 500 m map.
#
# Parte B - leave-station-out on N_DOBRAS official stations:
#   ccd_frio (as in the article), ccd_quente (start from the full-data
#   mode), eb_quente, theta_fixo. Time per fold, RMSE, CRPS, coverage,
#   and how far each moves from ccd_frio.
#
# Parte C - prediction grid at RESOLUCOES (from the 01 fit, no refit):
#   time with and without posterior samples, and how the coarse map
#   compares with the 500 m map averaged up to the same cells.
#
# Run alone on the machine (other jobs distort the times), e.g. after
# 06_sensibilidade_malha.R. Each step is saved as it finishes; rerunning
# skips what is done. Rode de dentro de analises/rio.
#   Rscript 07_custo.R > log_custo.log 2>&1
# Time: Parte A ~2-3 h, Parte B ~4-6 h (N_DOBRAS = 6), Parte C < 1 h.
# ==================================================================

library(bruce)
library(terra)
library(sf)
suppressMessages(library(INLA))
source("config.R")
inla.setOption(fmesher.evolution.warn = FALSE)

VAR <- "median"; MES <- "2025-01"
RODAR_A <- F; RODAR_B <- F; RODAR_C <- TRUE
N_DOBRAS   <- 6
RESOLUCOES <- c(500, 1000, 2000, 5000)
N_AMOSTRAS <- 200
DIR <- file.path(RES, "custo")
dir.create(DIR, showWarnings = FALSE, recursive = TRUE)

novo <- readRDS(sprintf(file.path(RES, "bruce_%s_%s.rds"), VAR, MES))
fit0 <- novo$fit
d    <- fit0$data
mesh <- fit0$mesh
th0  <- fit0$inla$mode$theta          # mode of the 01 fit, for warm starts / fixed theta
spec0 <- bruce_spec(range_prior = 15000, intercept_prior = c(0, 0.001), on_collapse = "warn")
janela <- range(fit0$days)

mrj <- geobr::read_municipality(code_muni = 3304557, year = 2020) |> st_transform(31983)
grid500 <- bruce_grid(mrj, 500, crs = fit0)
rs500 <- c(elev = "average", urb = "near", svf = "near")

nc <- parallel::detectCores()
cat(sprintf("Maquina: %d nucleos; INLA %s; malha %d vertices; %d dias\n",
            nc, packageVersion("INLA"), mesh$n, length(fit0$days)))

minutos <- function(t0) (proc.time()[["elapsed"]] - t0) / 60
hyper <- function(fit, nome) {
  hp <- fit$inla$summary.hyperpar
  if (is.null(hp) || !nome %in% rownames(hp)) return(NA_real_)
  hp[nome, "mean"]
}

# ==================================================================
# Parte A - one fit, several ways to run INLA
# ==================================================================
if (RODAR_A) {
  pardiso_ok <- tryCatch(any(grepl("SUCCESS", inla.pardiso.check())), error = function(e) FALSE)
  if (!pardiso_ok) cat("PARDISO sem licenca nesta maquina: configuracao 'pardiso' pulada\n",
                       "  (licenca academica gratuita: https://panua.ch/pardiso/ ; depois\n",
                       "   inla.setOption(pardiso.license = '~/.pardiso.lic'))\n")

  com <- function(...) { s <- spec0; a <- list(...); s[names(a)] <- a; s }
  configs <- list(
    referencia  = list(spec = spec0),
    eb          = list(spec = com(int_strategy = "eb")),
    eb_gaussian = list(spec = com(int_strategy = "eb", strategy = "gaussian")),
    theta_fixo  = list(spec = spec0, theta_fixo = TRUE)
  )
  for (tt in unique(c(sprintf("%d:1", nc), sprintf("%d:2", max(1, nc %/% 2))))) {
    configs[[paste0("threads_", sub(":", "x", tt))]] <- list(spec = com(num_threads = tt))
  }
  if (pardiso_ok) configs$pardiso <- list(spec = spec0, pardiso = TRUE)

  smtp0 <- inla.getOption("smtp")
  for (nm in names(configs)) {
    arq <- file.path(DIR, sprintf("A_%s.rds", nm))
    if (file.exists(arq)) { cat(sprintf("[A %s] ja existe\n", nm)); next }
    cf <- configs[[nm]]
    cat(sprintf("\n===== A: %s =====\n", nm))
    if (isTRUE(cf$pardiso)) inla.setOption(smtp = "pardiso")
    t0 <- proc.time()[["elapsed"]]
    fit <- tryCatch(
      if (isTRUE(cf$theta_fixo)) bruce_fit(d, mesh, cf$spec, window = janela, theta_init = th0, fix_theta = TRUE)
      else bruce_fit(d, mesh, cf$spec, window = janela),
      error = function(e) { message("ERRO: ", conditionMessage(e)); NULL })
    m <- minutos(t0)
    inla.setOption(smtp = smtp0)
    if (is.null(fit)) next
    pred <- bruce_predict(fit, grid = grid500, cap = FALSE, resample = rs500)
    saveRDS(list(nome = nm, minutos = m, cpu = fit$inla$cpu.used, status = fit$status,
                 range = hyper(fit, "Range for s"), stdev = hyper(fit, "Stdev for s"),
                 grouprho = hyper(fit, "GroupRho for s"),
                 fixed = fit$inla$summary.fixed, sd_campo = mean(fit$inla$summary.random$s$sd),
                 pred_mean = wrap(pred$mean)), arq)
    cat(sprintf("[A %s] %.1f min, status %s\n", nm, m, fit$status))
    rm(fit, pred); gc()
  }

  res <- lapply(list.files(DIR, "^A_.*\\.rds$", full.names = TRUE), readRDS)
  names(res) <- vapply(res, `[[`, "", "nome")
  if (!is.null(res$referencia)) {
    ref <- res$referencia
    mref <- unwrap(ref$pred_mean)
    tabA <- do.call(rbind, lapply(res, function(r) {
      fx <- r$fixed[rownames(ref$fixed), ]
      dm <- values(unwrap(r$pred_mean) - mref, mat = FALSE)
      data.frame(config = r$nome, minutos = round(r$minutos, 1),
                 acelera = round(ref$minutos / r$minutos, 1), status = r$status,
                 range = round(r$range), stdev = round(r$stdev, 2), grouprho = round(r$grouprho, 2),
                 max_dif_fixos_em_sd = round(max(abs(fx$mean - ref$fixed$mean) / ref$fixed$sd), 3),
                 razao_sd_campo = round(r$sd_campo / ref$sd_campo, 3),
                 rmse_mapa = round(sqrt(mean(dm^2, na.rm = TRUE)), 4),
                 max_dif_mapa = round(max(abs(dm), na.rm = TRUE), 3))
    }))
    tabA <- tabA[order(-tabA$acelera), ]; rownames(tabA) <- NULL
    cat("\n=== Parte A: custo do ajuste (mapas de 500 m comparados com a referencia) ===\n")
    print(tabA)
    write.csv(tabA, file.path(DIR, "resumo_A_ajuste.csv"), row.names = FALSE)
  }
}

# ==================================================================
# Parte B - leave-station-out on a subset of official stations
# ==================================================================
if (RODAR_B) {
  ofic <- d$obs[d$obs$source == "official" & d$obs$date %in% fit0$days, ]
  n_dias <- table(ofic$station_id)
  candidatas <- sort(names(n_dias)[n_dias >= 20])
  set.seed(42)
  ids <- sort(sample(candidatas, min(N_DOBRAS, length(candidatas))))
  cat("\nParte B: estacoes de teste:", paste(ids, collapse = ", "), "\n")

  configsB <- list(
    ccd_frio   = list(spec = spec0, theta = NULL, fixo = FALSE),
    ccd_quente = list(spec = spec0, theta = th0,  fixo = FALSE),
    eb_quente  = list(spec = { s <- spec0; s$int_strategy <- "eb"; s }, theta = th0, fixo = FALSE),
    theta_fixo = list(spec = spec0, theta = th0,  fixo = TRUE)
  )
  cvs <- list()
  for (nm in names(configsB)) {
    cf <- configsB[[nm]]
    cat(sprintf("\n===== B: %s =====\n", nm))
    cvs[[nm]] <- tryCatch(
      bruce_cv(d, mesh, cf$spec, window = janela, scheme = "lso", fold_ids = ids,
               n_samples = N_AMOSTRAS, seed = 42, theta_init = cf$theta, fix_theta = cf$fixo,
               checkpoint_dir = file.path(DIR, paste0("B_", nm))),
      error = function(e) { message("ERRO: ", conditionMessage(e)); NULL })
  }
  cvs <- Filter(Negate(is.null), cvs)
  if (!is.null(cvs$ccd_frio)) {
    base <- as.data.frame(cvs$ccd_frio$predictions)
    tabB <- do.call(rbind, lapply(names(cvs), function(nm) {
      cv <- cvs[[nm]]; o <- as.data.frame(cv$overall)
      p <- merge(base[, c("station_id", "date", "pred", "pred_sd")],
                 as.data.frame(cv$predictions)[, c("station_id", "date", "pred", "pred_sd")],
                 by = c("station_id", "date"), suffixes = c("_ref", ""))
      data.frame(config = nm, min_por_dobra = round(mean(cv$folds$minutes, na.rm = TRUE), 1),
                 rmse = round(o$rmse, 3), crps = round(o$crps, 3), cobertura90 = round(o$coverage, 3),
                 dif_media_pred = round(mean(abs(p$pred - p$pred_ref)), 3),
                 razao_sd_pred = round(mean(p$pred_sd) / mean(p$pred_sd_ref), 3),
                 status_nao_ok = sum(cv$folds$status != "ok"))
    }))
    cat("\n=== Parte B: LSO em", length(ids), "estacoes ===\n"); print(tabB)
    write.csv(tabB, file.path(DIR, "resumo_B_cv.csv"), row.names = FALSE)
  }
}

# ==================================================================
# Parte C - coarser prediction grids (01 fit, no refit)
# ==================================================================
if (RODAR_C) {
  rsC <- c(elev = "average", urb = "average", svf = "average")   # cell means of the covariates
  for (res_m in RESOLUCOES) {
    arq <- file.path(DIR, sprintf("C_%dm.rds", res_m))
    if (file.exists(arq)) { cat(sprintf("[C %d m] ja existe\n", res_m)); next }
    g <- bruce_grid(mrj, res_m, crs = fit0)
    rs <- if (res_m == 500) rs500 else rsC        # 500 m exactly as in 01
    t0 <- proc.time()[["elapsed"]]
    pm <- bruce_predict(fit0, grid = g, cap = FALSE, resample = rs)
    m_media <- minutos(t0)
    t0 <- proc.time()[["elapsed"]]
    ps <- bruce_predict(fit0, grid = g, cap = FALSE, resample = rs, n_samples = N_AMOSTRAS, seed = 42)
    m_amostras <- minutos(t0)
    saveRDS(list(res = res_m, celulas = length(pm$cells), min_media = m_media,
                 min_amostras = m_amostras, pred_mean = wrap(pm$mean)), arq)
    cat(sprintf("[C %d m] %d celulas; media %.2f min, %d amostras %.1f min\n",
                res_m, length(pm$cells), m_media, N_AMOSTRAS, m_amostras))
    rm(pm, ps); gc()
  }

  resC <- lapply(list.files(DIR, "^C_[0-9]+m\\.rds$", full.names = TRUE), readRDS)
  resC <- resC[order(vapply(resC, `[[`, 0, "res"))]
  m500 <- unwrap(resC[[which(vapply(resC, `[[`, 0, "res") == 500)]]$pred_mean)
  # Two ways to get a coarse map: predict directly on the coarse grid, or
  # predict at 500 m and average up. They differ because the spatial field
  # is evaluated at the cell centre in the first.
  tabC <- do.call(rbind, lapply(resC, function(r) {
    m <- unwrap(r$pred_mean)
    agreg <- resample(m500, m, method = "average")
    dif <- values(m - agreg, mat = FALSE)
    v <- values(m, mat = TRUE)
    data.frame(resolucao_m = r$res, celulas = r$celulas,
               min_media = round(r$min_media, 2), min_amostras = round(r$min_amostras, 1),
               sd_espacial = round(mean(apply(v, 2, sd, na.rm = TRUE)), 3),
               rmse_vs_500_agregado = round(sqrt(mean(dif^2, na.rm = TRUE)), 3),
               max_dif_vs_500_agregado = round(max(abs(dif), na.rm = TRUE), 2))
  }))
  cat("\n=== Parte C: resolucao da grade de predicao ===\n"); print(tabC)
  write.csv(tabC, file.path(DIR, "resumo_C_resolucao.csv"), row.names = FALSE)

  dia <- names(m500)[which.max(global(m500, "mean", na.rm = TRUE)$mean)]
  mapas <- lapply(resC, function(r) unwrap(r$pred_mean)[[dia]])
  lim <- range(unlist(lapply(mapas, function(m) range(values(m), na.rm = TRUE))))
  png(file.path(DIR, "resolucao.png"), width = 420 * length(mapas), height = 380, res = 110)
  par(mfrow = c(1, length(mapas)), mar = c(1, 1, 2.5, 3))
  for (k in seq_along(mapas)) {
    plot(mapas[[k]], col = hcl.colors(60, "YlOrRd", rev = TRUE), range = lim, axes = FALSE,
         main = sprintf("%g km", resC[[k]]$res / 1000))
    lines(vect(mrj), col = "grey40", lwd = 0.5)
  }
  dev.off()
  cat("Figura:", file.path(DIR, "resolucao.png"), "(dia", dia, ")\n")
}
