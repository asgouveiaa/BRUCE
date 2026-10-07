# ==================================================================
# 09_previsao_dia_seguinte.R
# Can bruce forecast the day after a fitted window, using only what it
# already knows (persistence through its day-to-day AR1 terms)?
#
# For each of the article's windows (e.g. January 2025), the model is fitted
# with the window extended by one day that has no data at all (1 February),
# with the window's own hyperparameters fixed (from 08_theta_emprestado.R),
# so each fit takes about a minute and never sees the day it forecasts. The
# forecast for that day is then compared, at the official stations and
# airports, with what they actually measured, and with two naive forecasts:
#   hoje        tomorrow = the station's value on the last day of the window
#   media_mes   tomorrow = the station's mean over the window
# The forecast is on the scale of a station reading: field + the station's
# own estimated offset + its network offset, with measurement noise in the
# interval.
#
# Also split into the city level (mean over stations) and the spatial
# pattern (deviation of each station from that mean): persistence is
# expected to hold better for the pattern than for the level.
#
# Needs: 08_theta_emprestado.R parte 1 done (resultados/theta_emprestado),
# the article packages (PACOTES) and query_meteorj_periodo.R for the
# observations of the forecast days (official stations and airports only;
# the Wunderground API is not called).
# Rode de dentro de analises/rio:
#   Rscript 09_previsao_dia_seguinte.R > log_previsao.log 2>&1
# Time: ~1-3 min per window (9 windows).
# ==================================================================

Sys.setenv(WU_API_KEY = "")          # cache only: the forecast check needs no Wunderground
library(bruce)
library(sf)
suppressMessages(library(INLA))
source("config.R")
source("query_meteorj_periodo.R")
inla.setOption(fmesher.evolution.warn = FALSE)
if (!"fix_theta" %in% names(formals(bruce_fit))) {
  stop("O bruce instalado e antigo (sem fix_theta). Reinstale: install.packages('~/bruce', repos = NULL, type = 'source')")
}

VAR <- "median"
PACOTES <- path.expand("~/MEGA/CIE things/pacotes_artigo")
DIR_THETA <- file.path(RES, "theta_emprestado")
DIR <- file.path(RES, "previsao_dia_seguinte")
dir.create(DIR, showWarnings = FALSE, recursive = TRUE)
N_AMOSTRAS <- 100
NIVEL <- 0.9
AVALIAR <- c("official", "airport")

arqs <- list.files(PACOTES, sprintf("^ARTIGO_pacote_%s_[0-9]{4}-[0-9]{2}\\.rds$", VAR), full.names = TRUE)
janelas <- sub(".*_([0-9]{4}-[0-9]{2})\\.rds$", "\\1", basename(arqs))
names(arqs) <- janelas
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

# The article's inputs for one window (same as 04 and 08)
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

# Forecast of a station READING on day g: field (samples) + network offset +
# the station's own estimated offset; sd adds that offset's sd and the
# network's measurement noise.
previsao_estacoes <- function(fit, st, g, hyper) {
  d <- fit$data
  cz <- d$covariates$station[match(st$station_id, d$covariates$station$station_id), fit$cov_names, drop = FALSE]
  X <- cbind(intercept = 1, as.matrix(cz))
  A <- INLA::inla.spde.make.A(fit$mesh, loc = cbind(st$x, st$y))
  comps <- bruce:::latent_components(fit, N_AMOSTRAS, NULL)
  P <- bruce:::eval_predictor(fit, comps, X, A, g)
  fx <- fit$inla$summary.fixed
  delta <- vapply(st$source, function(s) {
    col <- fit$bias_cols[s]; if (is.na(col)) 0 else fx[col, "mean"] }, 0)
  re <- fit$inla$summary.random$station_idx
  k <- match(st$station_id, fit$station_ids)
  sd_pop <- 1 / sqrt(hyper["Precision for station_idx", "mean"])   # from the window's own fit
  u_mean <- ifelse(is.na(k), 0, re$mean[k])
  u_sd <- ifelse(is.na(k), sd_pop, re$sd[k])
  ruido <- d$sources$obs_sd[match(st$source, d$sources$source)]
  data.frame(station_id = st$station_id, source = st$source,
             prev = rowMeans(P) + delta + u_mean,
             prev_sd = sqrt(apply(P, 1, stats::sd)^2 + u_sd^2 + ruido^2))
}

linhas <- list()
for (j in janelas) {
  arq <- file.path(DIR, sprintf("previsao_%s.rds", j))
  if (file.exists(arq)) { linhas[[j]] <- readRDS(arq); cat(sprintf("[%s] ja existe\n", j)); next }
  f_theta <- file.path(DIR_THETA, sprintf("proprio_%s.rds", j))
  if (!file.exists(f_theta)) { cat(sprintf("[%s] sem theta (rode a parte 1 do 08), pulado\n", j)); next }
  proprio <- readRDS(f_theta); th <- proprio$theta

  orig <- readRDS(arqs[[j]])
  d <- entradas(orig)
  ultimo <- max(d$obs$date); alvo <- ultimo + 1
  cat(sprintf("\n===== %s: previsao de %s =====\n", j, format(alvo)))

  # What the stations measured on the forecast day (official + airports)
  ob <- tryCatch(buscar_meteo_periodo_bulk(alvo, alvo, estatisticas = unique(c("median", "max", VAR))),
                 error = function(e) { message("ERRO na query: ", conditionMessage(e)); NULL })
  if (is.null(ob) || !nrow(ob)) { cat("  sem observacoes do dia previsto, pulado\n"); next }
  ob <- data.frame(station_id = ob$estacao_id, date = as.Date(ob$data),
                   obs = ob[[paste0(VAR, "_temperatura")]])
  ob <- ob[ob$date == alvo & !is.na(ob$obs), ]

  st <- unique(d$obs[d$obs$source %in% AVALIAR, c("station_id", "source", "x", "y")])
  st <- st[!duplicated(st$station_id) & st$station_id %in% ob$station_id, ]
  if (!nrow(st)) { cat("  nenhuma estacao oficial com dado no dia previsto, pulado\n"); next }

  t0 <- proc.time()[["elapsed"]]
  fit <- tryCatch(bruce_fit(d, orig$mesh, spec_janela(j), window = c(min(d$obs$date), alvo),
                            theta_init = th, fix_theta = TRUE, keep_config = TRUE),
                  error = function(e) { message("ERRO no ajuste: ", conditionMessage(e)); NULL })
  if (is.null(fit)) next
  pv <- previsao_estacoes(fit, st, which(fit$days == alvo), proprio$hyper)
  minutos <- (proc.time()[["elapsed"]] - t0) / 60

  hoje <- d$obs[d$obs$date == ultimo, c("station_id", "value")]
  media <- aggregate(value ~ station_id, d$obs, mean)
  res <- merge(pv, ob[, c("station_id", "obs")], by = "station_id")
  res$hoje <- hoje$value[match(res$station_id, hoje$station_id)]
  res$media_mes <- media$value[match(res$station_id, media$station_id)]
  res$janela <- j; res$dia <- alvo
  saveRDS(res, arq)
  linhas[[j]] <- res
  cat(sprintf("  %d estacoes, %.1f min | RMSE bruce %.2f, hoje %.2f, media do mes %.2f C\n",
              nrow(res), minutos, sqrt(mean((res$prev - res$obs)^2)),
              sqrt(mean((res$hoje - res$obs)^2, na.rm = TRUE)),
              sqrt(mean((res$media_mes - res$obs)^2, na.rm = TRUE))))
  rm(fit, orig, d); gc()
}

todos <- do.call(rbind, linhas)
if (is.null(todos) || !nrow(todos)) stop("Nenhuma previsao feita.")
write.csv(todos, file.path(DIR, "previsao_por_estacao.csv"), row.names = FALSE)
z <- qnorm(1 - (1 - NIVEL) / 2)

# Errors overall, split into city level and spatial pattern
avalia <- function(p, o) {
  ok <- !is.na(p) & !is.na(o)
  p <- p[ok]; o <- o[ok]
  c(rmse = sqrt(mean((p - o)^2)), vies = mean(p - o),
    erro_nivel = mean(p) - mean(o),
    rmse_padrao = sqrt(mean(((p - mean(p)) - (o - mean(o)))^2)),
    cor_padrao = if (length(o) > 2) stats::cor(p, o) else NA)
}
por_janela <- do.call(rbind, lapply(split(todos, todos$janela), function(r) {
  out <- lapply(c(bruce = "prev", hoje = "hoje", media_mes = "media_mes"),
                function(m) avalia(r[[m]], r$obs))
  data.frame(janela = r$janela[1], dia = format(r$dia[1]), n = nrow(r),
             rmse_bruce = out$bruce["rmse"], rmse_hoje = out$hoje["rmse"],
             rmse_media = out$media_mes["rmse"],
             nivel_bruce = out$bruce["erro_nivel"], nivel_hoje = out$hoje["erro_nivel"],
             padrao_bruce = out$bruce["rmse_padrao"], padrao_hoje = out$hoje["rmse_padrao"],
             cobertura = mean(abs(r$obs - r$prev) <= z * r$prev_sd))
}))
por_janela[, -(1:3)] <- round(por_janela[, -(1:3)], 2)
rownames(por_janela) <- NULL
cat("\n=== Previsao do dia seguinte, por janela (C) ===\n")
print(por_janela, row.names = FALSE)
cat("nivel_*: erro na media da cidade (previsto - observado)\n")
cat("padrao_*: RMSE depois de tirar a media da cidade (so o padrao espacial)\n")
cat(sprintf("cobertura: fracao das leituras dentro do intervalo de %.0f%% do bruce\n", 100 * NIVEL))
write.csv(por_janela, file.path(DIR, "previsao_por_janela.csv"), row.names = FALSE)

geral <- data.frame(
  metodo = c("bruce", "hoje (persistencia)", "media do mes"),
  rmse = c(sqrt(mean((todos$prev - todos$obs)^2)),
           sqrt(mean((todos$hoje - todos$obs)^2, na.rm = TRUE)),
           sqrt(mean((todos$media_mes - todos$obs)^2, na.rm = TRUE))),
  rmse_nivel = c(sqrt(mean(por_janela$nivel_bruce^2)), sqrt(mean(por_janela$nivel_hoje^2)),
                 sqrt(mean(vapply(split(todos, todos$janela), function(r)
                   mean(r$media_mes - r$obs, na.rm = TRUE), 0)^2))),
  rmse_padrao = c(sqrt(mean(por_janela$padrao_bruce^2)), sqrt(mean(por_janela$padrao_hoje^2)), NA))
geral[, -1] <- round(geral[, -1], 2)
cat(sprintf("\n=== Todas as janelas juntas (%d leituras, %d dias) ===\n", nrow(todos),
            length(unique(todos$janela))))
print(geral, row.names = FALSE)
cat(sprintf("Cobertura do intervalo de %.0f%% do bruce: %.2f\n", 100 * NIVEL,
            mean(abs(todos$obs - todos$prev) <= z * todos$prev_sd)))
write.csv(geral, file.path(DIR, "previsao_geral.csv"), row.names = FALSE)
