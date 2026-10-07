# ==================================================================
# 08_theta_emprestado.R
# Can the hyperparameters (theta) learned in one window be reused in
# another, and do the article's findings survive?
#
# Uses the article's own inputs for each of its windows (the final packages
# ARTIGO_pacote_<VAR>_<AAAA-MM>.rds in PACOTES: stations, ERA5 cells,
# standardised covariates and mesh), as in 04_comparar_entradas_identicas.R.
# Same model as the article in every window: where the article fixed the
# day-to-day field correlation (GroupRho), read from MANIFESTO_pacotes_artigo.csv,
# it is fixed here at the same value; if the article fixed the daily term's
# sd, give it in SD_DIARIO_FIXO (the script stops if it detects that case).
#
# Parte 1 - each window fitted by bruce with int_strategy = "eb", learning
#           its own theta. Compared with the article's fit (fixed effects,
#           hyperparameters): does "eb" reproduce the article?
# Parte 2 - each window fitted again with theta BORROWED from every other
#           window (fix_theta = TRUE). Compared with the window's own fit:
#           is elevation still significant? how much do the fixed effects
#           and the surface at the stations move?
#   Windows with different structure (e.g. a network present in one and
#   absent in the other) cannot share theta and are skipped.
#
# Rode de dentro de analises/rio (le RES de config.R). Each fit is
# saved as it finishes; rerunning skips what is done.
#   Rscript 08_theta_emprestado.R > log_theta.log 2>&1
# Time: Parte 1 = 9 fits with eb (~12 min each on the server; more on a
# laptop); Parte 2 = up to 72 fits with fixed theta (~1 min each).
# ==================================================================

library(bruce)
library(sf)
suppressMessages(library(INLA))
`%||%` <- function(a, b) if (is.null(a)) b else a
source("config.R")
if (!"fix_theta" %in% names(formals(bruce_fit))) {
  stop("O bruce instalado e antigo (sem fix_theta). Reinstale a partir de ~/bruce:\n",
       "  install.packages('~/bruce', repos = NULL, type = 'source')  e reinicie o R")
}
inla.setOption(fmesher.evolution.warn = FALSE)

VAR <- "median"
JANELAS <- NULL     # NULL = all windows of the article; or e.g. c("2024-11", "2025-01", "2025-06")
PACOTES <- path.expand("~/MEGA/CIE things/pacotes_artigo")   # final packages of the article
SD_DIARIO_FIXO <- NULL   # sd of the daily term, only if the article fixed it
DIR <- file.path(RES, "theta_emprestado")
dir.create(DIR, showWarnings = FALSE, recursive = TRUE)

arqs <- list.files(PACOTES, sprintf("^ARTIGO_pacote_%s_[0-9]{4}-[0-9]{2}\\.rds$", VAR),
                   full.names = TRUE)
if (!length(arqs)) stop("Nenhum pacote do artigo em ", PACOTES, " (confira PACOTES no topo)")
manifesto <- read.csv(file.path(PACOTES, "MANIFESTO_pacotes_artigo.csv"))
manifesto <- manifesto[manifesto$stat == VAR, ]
janelas <- sub(".*_([0-9]{4}-[0-9]{2})\\.rds$", "\\1", basename(arqs))
names(arqs) <- janelas
if (!is.null(JANELAS)) { janelas <- intersect(janelas, JANELAS); arqs <- arqs[janelas] }
cat("Janelas do artigo encontradas:", paste(janelas, collapse = ", "), "\n")

spec <- bruce_spec(range_prior = 15000, intercept_prior = c(0, 0.001), on_collapse = "warn",
                   int_strategy = "eb", day_sd_fixed = SD_DIARIO_FIXO)

# The article's model for window j: GroupRho fixed where the manifest says so
spec_janela <- function(j) {
  s <- spec
  m <- manifesto[manifesto$mes == j, ]
  if (nrow(m) == 1 && m$grouprho == "fixo") s$group_rho_fixed <- m$grouprho_valor
  s
}
grouprho_janela <- function(j) {
  m <- manifesto[manifesto$mes == j, ]
  if (nrow(m) == 1 && m$grouprho == "fixo") sprintf("fixo %.2f", m$grouprho_valor) else "livre"
}
# Stop if the article estimated no precision for the daily term (= it was fixed)
# and SD_DIARIO_FIXO was not given.
confere_diario <- function(orig, j) {
  estimado <- any(grepl("^Precision for dia_idx", rownames(orig$fit$summary.hyperpar)))
  if (!estimado && is.null(SD_DIARIO_FIXO)) {
    stop("Na janela ", j, " o artigo nao estimou a precisao do termo diario (dia_idx): ",
         "ela foi fixada. Informe o sd usado em SD_DIARIO_FIXO no topo do script.")
  }
  if (estimado && !is.null(SD_DIARIO_FIXO)) {
    stop("Na janela ", j, " o artigo estimou o termo diario; deixe SD_DIARIO_FIXO <- NULL.")
  }
}

# ---- The article's inputs for one window, on the bruce contract (as in 04) ----
entradas <- function(orig) {
  est <- st_drop_geometry(orig$est_all)
  src <- unname(c(ofic = "official", aura = "aurassure", wunder = "wunder")[est$origem])
  src[est$estacao_id %in% orig$ESTACOES_AEROPORTO] <- "airport"
  stations <- data.frame(station_id = est$estacao_id, source = src, x = est$x, y = est$y)
  stations <- stations[!duplicated(stations$station_id), ]
  obs <- data.frame(station_id = est$estacao_id, date = as.Date(est$data),
                    value = est[[orig$col_estacao]])
  bg <- st_drop_geometry(orig$d_all)
  background <- data.frame(x = bg$x, y = bg$y, date = as.Date(bg$date), value = bg[[orig$var_era5]])
  sources <- data.frame(source = c("official", "airport", "aurassure", "wunder"),
                        bias = c(FALSE, TRUE, TRUE, TRUE), obs_sd = c(0.3, 0.3, NA, NA))
  sources <- sources[sources$source %in% stations$source, ]
  d <- bruce_data(obs, stations, sources, crs = 31983, background = background, background_sd = 0.7)
  st_cov <- data.frame(station_id = est$estacao_id, elev = est$elev_std, urb = est$urb_std,
                       svf = est$svf_std)
  st_cov <- st_cov[!duplicated(st_cov$station_id), ]
  st_cov$svf[is.na(st_cov$svf)] <- 0            # INLA treats a missing covariate as 0
  bg_cov <- data.frame(x = bg$x, y = bg$y, elev = bg$elev_std, urb = bg$urb_std, svf = 0)
  bg_cov <- bg_cov[!duplicated(bg_cov[, c("x", "y")]), ]
  d$covariates <- list(names = c("elev", "urb", "svf"), rasters = NULL, scaling = NULL,
                       station = st_cov, background = bg_cov)
  d
}

# Posterior mean surface (reference scale) at every station location and day:
# stations x days. Used to measure how much borrowing theta moves the map.
superficie <- function(fit) {
  d <- fit$data
  st <- unique(d$obs[, c("station_id", "x", "y")])
  st <- st[!duplicated(st$station_id), ]
  cz <- d$covariates$station[match(st$station_id, d$covariates$station$station_id), fit$cov_names, drop = FALSE]
  X <- cbind(intercept = 1, as.matrix(cz))
  A <- INLA::inla.spde.make.A(fit$mesh, loc = cbind(st$x, st$y))
  comps <- bruce:::latent_components(fit, 0, NULL)
  m <- vapply(seq_along(fit$days), function(g) as.vector(bruce:::eval_predictor(fit, comps, X, A, g)),
              numeric(nrow(st)))
  dimnames(m) <- list(st$station_id, as.character(fit$days))
  m
}

fixos <- function(fx) fx[, c("mean", "sd", "0.025quant", "0.975quant")]
hp <- function(h, p) { i <- grep(p, rownames(h)); if (length(i)) h[i[1], "mean"] else NA_real_ }

# ==================================================================
# Parte 1 - own theta for every window (eb)
# ==================================================================
for (j in janelas) {
  arq <- file.path(DIR, sprintf("proprio_%s.rds", j))
  if (file.exists(arq)) { cat(sprintf("[proprio %s] ja existe\n", j)); next }
  cat(sprintf("\n===== Parte 1: %s =====\n", j))
  orig <- readRDS(arqs[[j]])
  confere_diario(orig, j)
  d <- entradas(orig)
  fit <- tryCatch(bruce_fit(d, orig$mesh, spec_janela(j)),
                  error = function(e) { message("ERRO: ", conditionMessage(e)); NULL })
  if (is.null(fit)) next
  saveRDS(list(janela = j, theta = fit$inla$mode$theta, status = fit$status, minutos = fit$minutes,
               fixed = fixos(fit$inla$summary.fixed), hyper = fit$inla$summary.hyperpar,
               sup = superficie(fit), fixed_artigo = fixos(orig$fit$summary.fixed),
               hyper_artigo = orig$fit$summary.hyperpar), arq)
  cat(sprintf("[proprio %s] %.1f min, status %s, range %.0f m\n", j, fit$minutes %||% NA_real_, fit$status,
              hp(fit$inla$summary.hyperpar, "^Range for s")))
  rm(fit, orig, d); gc()
}

prop <- lapply(janelas, function(j) {
  f <- file.path(DIR, sprintf("proprio_%s.rds", j)); if (file.exists(f)) readRDS(f) else NULL })
names(prop) <- janelas
prop <- Filter(Negate(is.null), prop)

# Does eb reproduce the article? (fixed effects in units of the article's sd)
nomes <- c(intercept = "b0", elev = "elev_std", urb = "urb_std", svf = "svf_std",
           bias_aurassure = "aura", bias_wunder = "wunder", bias_airport = "aeroporto")
tab1 <- do.call(rbind, lapply(prop, function(p) {
  nm <- nomes[names(nomes) %in% rownames(p$fixed) & nomes %in% rownames(p$fixed_artigo)]
  dif <- (p$fixed[names(nm), "mean"] - p$fixed_artigo[nm, "mean"]) / p$fixed_artigo[nm, "sd"]
  ea <- p$fixed_artigo["elev_std", ]; eb <- p$fixed["elev", ]
  data.frame(janela = p$janela, grouprho = grouprho_janela(p$janela), status = p$status,
             minutos = if (is.null(p$minutos)) NA else round(p$minutos, 1),
             range_artigo = round(hp(p$hyper_artigo, "^Range for s")),
             range_eb = round(hp(p$hyper, "^Range for s")),
             elev_artigo = sprintf("%.2f [%.2f; %.2f]", ea$mean, ea$`0.025quant`, ea$`0.975quant`),
             elev_eb = sprintf("%.2f [%.2f; %.2f]", eb$mean, eb$`0.025quant`, eb$`0.975quant`),
             max_dif_fixos_em_sd = round(max(abs(dif)), 2))
}))
cat("\n=== Parte 1: bruce (eb) contra o artigo, com as mesmas entradas ===\n")
print(tab1, row.names = FALSE)
write.csv(tab1, file.path(DIR, "parte1_eb_vs_artigo.csv"), row.names = FALSE)

# ==================================================================
# Parte 2 - every window with theta borrowed from every other window
# ==================================================================
doadores <- Filter(function(p) p$status == "ok", prop)
if (length(doadores) < length(prop)) {
  cat("Nao emprestam theta (ajuste nao 'ok'):",
      paste(setdiff(names(prop), names(doadores)), collapse = ", "), "\n")
}
for (j in names(prop)) {
  arq <- file.path(DIR, sprintf("emprestado_%s.rds", j))
  if (file.exists(arq)) { cat(sprintf("[emprestado %s] ja existe\n", j)); next }
  cat(sprintf("\n===== Parte 2: %s =====\n", j))
  orig <- readRDS(arqs[[j]])
  d <- entradas(orig)
  alvo <- prop[[j]]
  linhas <- list()
  for (s in setdiff(names(doadores), j)) {
    th <- doadores[[s]]$theta
    if (length(th) != length(alvo$theta) || !identical(names(th), names(alvo$theta))) {
      cat(sprintf("  theta de %s: estrutura diferente, pulado\n", s)); next
    }
    fit <- tryCatch(bruce_fit(d, orig$mesh, spec_janela(j), theta_init = th, fix_theta = TRUE),
                    error = function(e) { message("  ERRO com theta de ", s, ": ", conditionMessage(e)); NULL })
    if (is.null(fit)) next
    fx <- fixos(fit$inla$summary.fixed)
    comum <- intersect(rownames(fx), rownames(alvo$fixed))
    sup <- superficie(fit)
    dif_sup <- sup - alvo$sup[rownames(sup), colnames(sup)]
    linhas[[s]] <- data.frame(
      alvo = j, doador = s, minutos = round(fit$minutes, 2),
      elev = fx["elev", "mean"], elev_q025 = fx["elev", "0.025quant"], elev_q975 = fx["elev", "0.975quant"],
      elev_proprio = alvo$fixed["elev", "mean"],
      max_dif_fixos_em_sd = max(abs(fx[comum, "mean"] - alvo$fixed[comum, "mean"]) / alvo$fixed[comum, "sd"]),
      rmse_superficie = sqrt(mean(dif_sup^2)), max_dif_superficie = max(abs(dif_sup)))
    cat(sprintf("  theta de %s: elev %.2f [%.2f; %.2f], superficie RMSE %.3f C\n", s,
                fx["elev", "mean"], fx["elev", "0.025quant"], fx["elev", "0.975quant"],
                linhas[[s]]$rmse_superficie))
    rm(fit); gc()
  }
  saveRDS(do.call(rbind, linhas), arq)
  rm(orig, d); gc()
}

emp <- do.call(rbind, lapply(names(prop), function(j) {
  f <- file.path(DIR, sprintf("emprestado_%s.rds", j)); if (file.exists(f)) readRDS(f) }))
if (!is.null(emp) && nrow(emp)) {
  emp$elev_significativo <- emp$elev_q025 > 0 | emp$elev_q975 < 0
  write.csv(emp, file.path(DIR, "parte2_todos_os_pares.csv"), row.names = FALSE)

  tab2 <- do.call(rbind, lapply(split(emp, emp$alvo), function(e) data.frame(
    janela = e$alvo[1], doadores = nrow(e),
    elev_proprio = round(e$elev_proprio[1], 2),
    elev_emprestado = sprintf("%.2f a %.2f", min(e$elev), max(e$elev)),
    elev_signif = sprintf("%d de %d", sum(e$elev_significativo), nrow(e)),
    max_dif_fixos_em_sd = round(max(e$max_dif_fixos_em_sd), 2),
    rmse_superficie_mediana = round(median(e$rmse_superficie), 3),
    rmse_superficie_pior = round(max(e$rmse_superficie), 3),
    pior_doador = e$doador[which.max(e$rmse_superficie)])))
  cat("\n=== Parte 2: cada janela com theta emprestado das outras ===\n")
  print(tab2, row.names = FALSE)
  cat("rmse_superficie: diferenca (C) da superficie nas estacoes em relacao ao theta proprio\n")
  write.csv(tab2, file.path(DIR, "parte2_resumo.csv"), row.names = FALSE)

  # Heatmap: surface RMSE, target (rows) x donor (columns)
  js <- names(prop)
  M <- matrix(NA_real_, length(js), length(js), dimnames = list(js, js))
  M[cbind(match(emp$alvo, js), match(emp$doador, js))] <- emp$rmse_superficie
  png(file.path(DIR, "rmse_superficie_alvo_x_doador.png"), width = 900, height = 800, res = 110)
  par(mar = c(5, 5, 3, 6))
  pal <- hcl.colors(30, "YlOrRd", rev = TRUE)
  image(seq_along(js), seq_along(js), t(M), col = pal, axes = FALSE, xlab = "theta de (doador)",
        ylab = "janela ajustada (alvo)", main = "Diferenca na superficie (RMSE, C)")
  axis(1, seq_along(js), js, las = 2, cex.axis = 0.8); axis(2, seq_along(js), js, las = 1, cex.axis = 0.8)
  ok <- which(!is.na(t(M)), arr.ind = TRUE)
  text(ok[, 1], ok[, 2], sprintf("%.2f", t(M)[ok]), cex = 0.65)
  dev.off()
  cat("Figura:", file.path(DIR, "rmse_superficie_alvo_x_doador.png"), "\n")
}

# How much does theta itself vary between windows?
tab3 <- do.call(rbind, lapply(prop, function(p) data.frame(
  janela = p$janela, grouprho_artigo = grouprho_janela(p$janela),
  range = round(hp(p$hyper, "^Range for s")),
  stdev_campo = round(hp(p$hyper, "^Stdev for s"), 2),
  grouprho = round(hp(p$hyper, "^GroupRho for s"), 2),
  rho_diario = round(hp(p$hyper, "^Rho for day_idx"), 2))))
cat("\n=== Hiperparametros de cada janela (eb) ===\n"); print(tab3, row.names = FALSE)
write.csv(tab3, file.path(DIR, "hiperparametros_por_janela.csv"), row.names = FALSE)
