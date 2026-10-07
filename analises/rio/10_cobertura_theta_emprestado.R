# ==================================================================
# 10_cobertura_theta_emprestado.R
# Are the uncertainty intervals still right when the hyperparameters
# (theta) are borrowed from another window? This is what a daily map
# built with last month's theta would rely on.
#
# Leave-station-out on N_DOBRAS official stations of a few target windows
# of the article (same inputs as 08), with:
#   proprio_eb     reference: every fold re-learns theta (eb, warm start
#                  from the window's own theta) - honest, no shortcut
#   proprio_fixo   theta of the window itself, fixed in every fold
#   de_<janela>    theta BORROWED from another window, fixed in every fold
# Compared on error (RMSE, CRPS) and on coverage of the 90% interval, both
# on the field scale and on the scale of a station reading (adding the
# station offset sd and the 0.3 C measurement noise, as in 07).
#
# Needs 08_theta_emprestado.R parte 1 (resultados/theta_emprestado).
# Rode de dentro de analises/rio:
#   Rscript 10_cobertura_theta_emprestado.R > log_cobertura.log 2>&1
# Time: ~45-60 min per target window with N_DOBRAS = 8 (the eb reference
# is most of it). Each fold is saved; rerunning resumes.
# ==================================================================

library(bruce)
library(sf)
suppressMessages(library(INLA))
source("config.R")

# These scripts score the field scale, as in the article (bruce >= 0.1.0.9
# scores station readings by default; older versions have no `scale`).
cv_campo <- function(...) {
  if ("scale" %in% names(formals(bruce_cv))) bruce_cv(..., scale = "field") else bruce_cv(...)
}
inla.setOption(fmesher.evolution.warn = FALSE)
if (!all(c("fix_theta", "fold_ids") %in% c(names(formals(bruce_fit)), names(formals(bruce_cv))))) {
  stop("O bruce instalado e antigo. Reinstale: install.packages('~/bruce', repos = NULL, type = 'source')")
}

VAR <- "median"
PACOTES <- path.expand("~/MEGA/CIE things/pacotes_artigo")
DIR_THETA <- file.path(RES, "theta_emprestado")
DIR <- file.path(RES, "cobertura_theta_emprestado")
dir.create(DIR, showWarnings = FALSE, recursive = TRUE)
N_DOBRAS <- 8
N_AMOSTRAS <- 100
NIVEL <- 0.9
# target window -> donors (pairs with the same structure, from 08). The last
# one is the worst pair seen in 08 (ranges 15 and 37 km).
ALVOS <- list(`2025-01` = c("2023-01", "2023-06", "2023-11"),
              `2024-11` = c("2024-01", "2024-06"),
              `2025-11` = c("2025-06"))

manifesto <- read.csv(file.path(PACOTES, "MANIFESTO_pacotes_artigo.csv"))
manifesto <- manifesto[manifesto$stat == VAR, ]
spec <- bruce_spec(range_prior = 15000, intercept_prior = c(0, 0.001), on_collapse = "warn",
                   int_strategy = "eb")
spec_janela <- function(j) {
  s <- spec
  m <- manifesto[manifesto$mes == j, ]
  if (nrow(m) == 1 && m$grouprho == "fixo") s$group_rho_fixed <- m$grouprho_valor
  s
}

# The article's inputs for one window (same as 04, 08 and 09)
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
  st_cov$svf[is.na(st_cov$svf)] <- 0
  bg_cov <- data.frame(x = bg$x, y = bg$y, elev = bg$elev_std, urb = bg$urb_std, svf = 0)
  bg_cov <- bg_cov[!duplicated(bg_cov[, c("x", "y")]), ]
  d$covariates <- list(names = c("elev", "urb", "svf"), rasters = NULL, scaling = NULL,
                       station = st_cov, background = bg_cov)
  d
}

theta_de <- function(j) {
  f <- file.path(DIR_THETA, sprintf("proprio_%s.rds", j))
  if (!file.exists(f)) stop("Falta ", f, " (rode a parte 1 do 08_theta_emprestado.R)")
  readRDS(f)
}

z <- qnorm(1 - (1 - NIVEL) / 2)
linhas <- list()
for (alvo in names(ALVOS)) {
  cat(sprintf("\n===== Janela %s =====\n", alvo))
  arq_pac <- file.path(PACOTES, sprintf("ARTIGO_pacote_%s_%s.rds", VAR, alvo))
  orig <- readRDS(arq_pac)
  d <- entradas(orig)
  proprio <- theta_de(alvo)
  sd_estacao <- 1 / sqrt(proprio$hyper["Precision for station_idx", "mean"])

  ofic <- d$obs[d$obs$source == "official", ]
  n_dias <- table(ofic$station_id)
  set.seed(42)
  cand <- sort(names(n_dias)[n_dias >= 20])
  ids <- sort(sample(cand, min(N_DOBRAS, length(cand))))
  cat("Estacoes de teste:", paste(ids, collapse = ", "), "\n")

  configs <- list(proprio_eb = list(theta = proprio$theta, fixo = FALSE),
                  proprio_fixo = list(theta = proprio$theta, fixo = TRUE))
  for (s in ALVOS[[alvo]]) {
    th <- theta_de(s)$theta
    if (length(th) != length(proprio$theta) || !identical(names(th), names(proprio$theta))) {
      cat(sprintf("  theta de %s: estrutura diferente, pulado\n", s)); next
    }
    configs[[paste0("de_", s)]] <- list(theta = th, fixo = TRUE)
  }

  cvs <- list()
  for (nm in names(configs)) {
    cat(sprintf("\n--- %s: %s ---\n", alvo, nm))
    cf <- configs[[nm]]
    cvs[[nm]] <- tryCatch(
      cv_campo(d, orig$mesh, spec_janela(alvo), scheme = "lso", fold_ids = ids,
               n_samples = N_AMOSTRAS, seed = NULL, theta_init = cf$theta, fix_theta = cf$fixo,
               checkpoint_dir = file.path(DIR, sprintf("%s_%s", alvo, nm))),
      error = function(e) { message("ERRO: ", conditionMessage(e)); NULL })
  }
  cvs <- Filter(Negate(is.null), cvs)
  if (is.null(cvs$proprio_eb)) next
  ref <- as.data.frame(cvs$proprio_eb$predictions)

  for (nm in names(cvs)) {
    p <- as.data.frame(cvs[[nm]]$predictions)
    m <- merge(ref[, c("station_id", "date", "pred", "pred_sd")],
               p[, c("station_id", "date", "pred", "pred_sd")],
               by = c("station_id", "date"), suffixes = c("_ref", ""))
    sd_obs <- sqrt(p$pred_sd^2 + sd_estacao^2 + 0.3^2)
    linhas[[paste(alvo, nm)]] <- data.frame(
      janela = alvo, config = nm, n = nrow(p),
      min_por_dobra = round(mean(cvs[[nm]]$folds$minutes, na.rm = TRUE), 1),
      rmse = round(sqrt(mean((p$pred - p$obs)^2)), 3),
      crps = round(mean(p$crps), 3),
      cob_campo = round(mean(p$obs >= p$lower & p$obs <= p$upper), 3),
      cob_leitura = round(mean(abs(p$obs - p$pred) <= z * sd_obs), 3),
      razao_sd = round(mean(m$pred_sd) / mean(m$pred_sd_ref), 3),
      dif_media_pred = round(mean(abs(m$pred - m$pred_ref)), 3))
  }
  rm(orig, d); gc()
}

tab <- do.call(rbind, linhas); rownames(tab) <- NULL
cat("\n=== Cobertura e erro com theta proprio e emprestado (LSO) ===\n")
print(tab, row.names = FALSE)
cat(sprintf("cob_campo / cob_leitura: fracao dentro do intervalo de %.0f%% (ideal ~%.2f na escala da leitura)\n",
            100 * NIVEL, NIVEL))
cat("razao_sd e dif_media_pred: em relacao a proprio_eb (que reaprende theta em cada dobra)\n")
write.csv(tab, file.path(DIR, "cobertura_theta_emprestado.csv"), row.names = FALSE)
