# ==================================================================
# 18_temperatura_comparar_simples.R
# Is bruce worth it for TEMPERATURE? Same question and same simple methods
# as 17 (dew point), now on the article's windows:
#
#   era5            ERA5-Land of the nearest cell, as it is
#   era5_vies       + the day's mean difference at the other official stations
#   era5_cov        ERA5 + regression on all the other stations: day + network
#                   + covariate differences (elevation, urban, sky view
#                   factor; station minus its ERA5 cell)
#   idw             inverse-distance weighting of the other stations'
#                   readings, citizen networks corrected by their offset
#   era5_cov_idw    era5_cov + IDW of its residuals (simple "covariates +
#                   spatial field")
#   bruce           leave-station-out with bruce, same stations and days
#
# Held out: the official stations (reference network) with >= 20 days, one
# at a time (up to N_DOBRAS per window). bruce re-learns theta in every fold
# starting from the window's own theta (eb, warm start: honest, as the
# reference of script 10); rapido=1 keeps that theta fixed instead (much
# faster, same error in script 10).
#
# Needs the article packages and, for the warm start, 08 parte 1
# (resultados/theta_emprestado/proprio_<window>.rds; without it bruce starts
# from scratch, slower). Every fold is saved: rerunning resumes.
#
# Rode de dentro de analises/rio (all windows: run it overnight):
#   Rscript 18_temperatura_comparar_simples.R > log_temp_comparar.log 2>&1 &
#   Rscript 18_temperatura_comparar_simples.R janelas=2025-01,2023-11 rapido=1
# Time: eb ~5-8 min per fold (~1 h per window with 8 folds); rapido=1 ~1 min.
# ==================================================================

library(bruce)
library(sf)
suppressMessages(library(INLA))
source("config.R")
inla.setOption(fmesher.evolution.warn = FALSE)

args <- commandArgs(trailingOnly = TRUE)
arg <- function(nome, padrao) {
  a <- grep(paste0("^", nome, "="), args, value = TRUE)
  if (length(a)) sub("^[^=]*=", "", a[1]) else padrao
}
VAR <- arg("var", "median")
RAPIDO <- arg("rapido", "0") == "1"
N_DOBRAS <- as.integer(arg("dobras", 8))
N_AMOSTRAS <- 100
POT_IDW <- 2
REF <- "official"
PACOTES <- path.expand("~/MEGA/CIE things/pacotes_artigo")
DIR_THETA <- file.path(RES, "theta_emprestado")
DIR <- file.path(RES, "comparar_simples_temperatura")
dir.create(DIR, showWarnings = FALSE, recursive = TRUE)

arqs <- list.files(PACOTES, sprintf("^ARTIGO_pacote_%s_[0-9]{4}-[0-9]{2}\\.rds$", VAR), full.names = TRUE)
names(arqs) <- sub(".*_([0-9]{4}-[0-9]{2})\\.rds$", "\\1", basename(arqs))
janelas <- if (arg("janelas", "") == "") names(arqs) else strsplit(arg("janelas", ""), ",")[[1]]
janelas <- intersect(janelas, names(arqs))
cat("Variavel:", VAR, " janelas:", paste(janelas, collapse = ", "), " modo:", if (RAPIDO) "theta fixo" else "eb com partida quente", "\n")

manifesto <- read.csv(file.path(PACOTES, "MANIFESTO_pacotes_artigo.csv"))
manifesto <- manifesto[manifesto$stat == VAR, ]
spec <- bruce_spec(range_prior = 15000, intercept_prior = c(0, 0.001), on_collapse = "warn", int_strategy = "eb")
spec_janela <- function(j) {
  s <- spec
  m <- manifesto[manifesto$mes == j, ]
  if (nrow(m) == 1 && m$grouprho == "fixo") s$group_rho_fixed <- m$grouprho_valor
  s
}

# The article's inputs for one window (same as 04, 08, 09, 10)
entradas <- function(orig) {
  est <- st_drop_geometry(orig$est_all)
  src <- unname(c(ofic = "official", aura = "aurassure", wunder = "wunder")[est$origem])
  src[est$estacao_id %in% orig$ESTACOES_AEROPORTO] <- "airport"
  stations <- data.frame(station_id = est$estacao_id, source = src, x = est$x, y = est$y)
  stations <- stations[!duplicated(stations$station_id), ]
  obs <- data.frame(station_id = est$estacao_id, date = as.Date(est$data), value = est[[orig$col_estacao]])
  bg <- st_drop_geometry(orig$d_all)
  background <- data.frame(x = bg$x, y = bg$y, date = as.Date(bg$date), value = bg[[orig$var_era5]])
  sources <- data.frame(source = c("official", "airport", "aurassure", "wunder"),
                        bias = c(FALSE, TRUE, TRUE, TRUE), obs_sd = c(0.3, 0.3, NA, NA))
  sources <- sources[sources$source %in% stations$source, ]
  d <- bruce_data(obs, stations, sources, crs = 31983, background = background, background_sd = 0.7)
  st_cov <- data.frame(station_id = est$estacao_id, elev = est$elev_std, urb = est$urb_std, svf = est$svf_std)
  st_cov <- st_cov[!duplicated(st_cov$station_id), ]
  st_cov$svf[is.na(st_cov$svf)] <- 0
  bg_cov <- data.frame(x = bg$x, y = bg$y, elev = bg$elev_std, urb = bg$urb_std, svf = 0)
  bg_cov <- bg_cov[!duplicated(bg_cov[, c("x", "y")]), ]
  d$covariates <- list(names = c("elev", "urb", "svf"), rasters = NULL, scaling = NULL,
                       station = st_cov, background = bg_cov)
  d
}

idw <- function(x0, y0, x, y, v) {
  w <- 1 / pmax(sqrt((x - x0)^2 + (y - y0)^2), 100)^POT_IDW
  sum(w * v) / sum(w)
}

# The simple methods, leaving each held-out station out, on the same days as bruce
simples <- function(d, pb) {
  obs <- as.data.frame(d$obs)[, c("station_id", "date", "value", "source", "x", "y")]
  bg <- as.data.frame(d$background)
  cel <- unique(bg[, c("x", "y")])
  bgc <- d$covariates$background
  cel <- merge(cel, bgc, by = c("x", "y"), all.x = TRUE, sort = FALSE)
  est <- unique(obs[, c("station_id", "x", "y")])
  est <- merge(est, d$covariates$station, by = "station_id", all.x = TRUE)
  est$cel <- vapply(seq_len(nrow(est)), function(i) which.min((cel$x - est$x[i])^2 + (cel$y - est$y[i])^2), 1L)
  est$d_elev <- est$elev - cel$elev[est$cel]
  est$d_urb <- est$urb - cel$urb[est$cel]
  est$d_svf <- est$svf
  bg$cel <- match(paste(bg$x, bg$y), paste(cel$x, cel$y))
  obs <- merge(obs, est[, c("station_id", "cel", "d_elev", "d_urb", "d_svf")], by = "station_id")
  obs$era5 <- bg$value[match(paste(obs$cel, obs$date), paste(bg$cel, bg$date))]
  obs <- obs[!is.na(obs$era5), ]
  obs$dif <- obs$value - obs$era5
  linhas <- list()
  for (s in unique(pb$station_id)) {
    tr <- obs[obs$station_id != s, ]
    te <- obs[obs$station_id == s & obs$date %in% pb$date[pb$station_id == s], ]
    if (!nrow(te)) next
    tr$rede <- relevel(factor(ifelse(tr$source == REF, "ref", tr$source)), "ref")
    m <- lm(dif ~ factor(date) + rede + d_elev + d_urb + d_svf, data = tr, na.action = na.exclude)
    off <- c(ref = 0, coef(m)[grep("^rede", names(coef(m)))])
    names(off) <- sub("^rede", "", names(off))
    off[is.na(off)] <- 0
    tr$corr <- tr$value - off[ifelse(tr$source == REF, "ref", tr$source)]
    tr$res <- residuals(m)            # NA where a covariate is missing (na.exclude)
    for (k in seq_len(nrow(te))) {
      dd <- te$date[k]
      ref_d <- tr[tr$date == dd & tr$source == REF, ]
      all_d <- tr[tr$date == dd & !is.na(tr$res), ]
      p_cov <- tryCatch(predict(m, newdata = data.frame(date = dd, rede = factor("ref", levels = levels(tr$rede)),
                                                        d_elev = te$d_elev[k], d_urb = te$d_urb[k], d_svf = te$d_svf[k])),
                        error = function(e) NA_real_)
      linhas[[length(linhas) + 1]] <- data.frame(
        station_id = s, date = dd, obs = te$value[k], era5 = te$era5[k],
        era5_vies = te$era5[k] + if (nrow(ref_d)) mean(ref_d$dif) else NA,
        era5_cov = te$era5[k] + p_cov,
        idw = if (nrow(all_d)) idw(te$x[k], te$y[k], all_d$x, all_d$y, all_d$corr) else NA,
        era5_cov_idw = te$era5[k] + p_cov + if (nrow(all_d)) idw(te$x[k], te$y[k], all_d$x, all_d$y, all_d$res) else 0)
    }
  }
  do.call(rbind, linhas)
}

metodos <- c("era5", "era5_vies", "era5_cov", "idw", "era5_cov_idw", "bruce")
nota <- function(p) do.call(rbind, lapply(metodos, function(mt) {
  e <- p[[mt]] - p$obs
  data.frame(metodo = mt, n = length(e), rmse = round(sqrt(mean(e^2)), 3),
             mae = round(mean(abs(e)), 3), vies = round(mean(e), 3))
}))

todas <- list()
for (j in janelas) {
  cat(sprintf("\n===== %s %s =====\n", VAR, j))
  arq_res <- file.path(DIR, sprintf("comparar_%s_%s%s.rds", VAR, j, if (RAPIDO) "_rapido" else ""))
  if (file.exists(arq_res)) { pr <- readRDS(arq_res); cat("ja existe, carregado\n") } else {
    orig <- readRDS(arqs[[j]])
    d <- entradas(orig)
    ofic <- d$obs[d$obs$source == REF, ]
    n_dias <- table(ofic$station_id)
    set.seed(42)
    cand <- sort(names(n_dias)[n_dias >= 20])
    ids <- sort(if (length(cand) > N_DOBRAS) sample(cand, N_DOBRAS) else cand)
    cat("Estacoes deixadas de fora:", paste(ids, collapse = ", "), "\n")
    f_theta <- file.path(DIR_THETA, sprintf("proprio_%s.rds", j))
    th <- if (VAR == "median" && file.exists(f_theta)) readRDS(f_theta)$theta else NULL
    if (is.null(th)) cat("(sem theta do 08: bruce comeca do zero em cada dobra)\n")
    cv_args <- list(d = d, mesh = orig$mesh, spec = spec_janela(j), scheme = "lso", fold_ids = ids,
                    n_samples = N_AMOSTRAS, seed = NULL, theta_init = th, fix_theta = RAPIDO && !is.null(th),
                    checkpoint_dir = file.path(DIR, sprintf("cv_%s_%s%s", VAR, j, if (RAPIDO) "_rapido" else "")))
    if ("scale" %in% names(formals(bruce_cv))) cv_args$scale <- "field"
    cv <- tryCatch(do.call(bruce_cv, cv_args), error = function(e) { message("ERRO bruce: ", conditionMessage(e)); NULL })
    if (is.null(cv)) next
    pb <- as.data.frame(cv$predictions)[, c("station_id", "date", "obs", "pred")]
    pb$date <- as.Date(pb$date)
    pr <- simples(d, pb)
    pr <- merge(pr, pb[, c("station_id", "date", "pred")], by = c("station_id", "date"))
    names(pr)[names(pr) == "pred"] <- "bruce"
    pr <- pr[stats::complete.cases(pr[, c("obs", metodos)]), ]
    pr$janela <- j
    saveRDS(pr, arq_res)
    rm(orig, d); gc()
  }
  g <- nota(pr)
  print(g[order(g$rmse), ], row.names = FALSE)
  todas[[j]] <- pr
}

if (length(todas)) {
  tudo <- do.call(rbind, todas)
  cat(sprintf("\n=== %s, todas as janelas juntas (%d janelas) ===\n", VAR, length(todas)))
  g <- nota(tudo)
  print(g[order(g$rmse), ], row.names = FALSE)

  cat("\nRMSE por janela:\n")
  pj <- do.call(rbind, lapply(split(tudo, tudo$janela), function(p) {
    v <- sapply(metodos, function(mt) round(sqrt(mean((p[[mt]] - p$obs)^2)), 2))
    data.frame(janela = p$janela[1], n = nrow(p), t(v))
  }))
  print(pj, row.names = FALSE)

  melhor <- g$metodo[g$metodo != "bruce"][which.min(g$rmse[g$metodo != "bruce"])]
  d_err <- abs(tudo$bruce - tudo$obs) - abs(tudo[[melhor]] - tudo$obs)
  # uncertainty by resampling whole stations (days of one station are not independent)
  set.seed(1)
  chave <- paste(tudo$janela, tudo$station_id)
  grupos <- split(d_err, chave)
  boot <- replicate(2000, mean(unlist(sample(grupos, replace = TRUE))))
  cat(sprintf("\nbruce x melhor metodo simples (%s): bruce erra menos em %.0f%% dos casos;\n",
              melhor, 100 * mean(d_err < 0)))
  cat(sprintf("diferenca media do erro absoluto %.3f C (IC95%% por estacao %.3f a %.3f; negativo = bruce melhor)\n",
              mean(d_err), quantile(boot, 0.025), quantile(boot, 0.975)))
  write.csv(tudo, file.path(DIR, sprintf("comparar_simples_temperatura_%s.csv", VAR)), row.names = FALSE)
  cat("\nSalvo em", file.path(DIR, sprintf("comparar_simples_temperatura_%s.csv", VAR)), "\n")
}
