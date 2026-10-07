# ==================================================================
# query_orvalho_periodo.R
# Daily dew point (Td) per station, for modelling with bruce.
#
# Same two-stage logic as the temperature query: hourly first (median of
# the records in the hour), then daily statistics over the hours.
#   INMET    Td measured (temperatura_orvalho); only the city's stations
#   REDEMET  Td measured (METAR, whole degrees); airports
#   Alerta Rio (tb_clima_ar)   Td computed from temperature + RH, record by record
#   Aurassure (tb_registers)   Td computed from temp + humid, record by record
#   Wunderground (API)         Td computed from tempAvg + humidityAvg, observation
#                              by observation; kept in its own cache
#                              (wu_orvalho_cache.duckdb, next to WU_CACHE), so the
#                              temperature cache is never touched. Needs
#                              WU_API_KEY for days not yet cached (~650 calls for
#                              21 stations x 31 days, once).
# Td is computed BEFORE any aggregation: the Magnus formula is non-linear, so
# Td from a daily median of T and RH is not the daily median of Td.
#
# Quality filters (each record): RH within [10, 100] %, T within [5, 45] C,
# Td within [0, 30] C and Td <= T + 0.5 C. The number of records removed is
# printed per source. For RH-based sources, `frac_saturada` gives the share of
# the day's hours with RH >= 99.5% (sensor stuck or true saturation; Td = T
# then): kept by default, filter with `max_frac_saturada` if needed.
# Days with fewer than `min_horas` hours are dropped (a median over a few
# hours is biased by the daily cycle).
#
# Use:
#   source("query_orvalho_periodo.R")
#   td <- buscar_orvalho_periodo("2026-01-01", "2026-01-31")
# ==================================================================

library(DBI)
library(duckdb)
library(dplyr)
library(glue)
library(httr2)
`%||%` <- function(a, b) if (is.null(a)) b else a

BASE_CLIMA <- "~/MEGA/CIE things/ducks/clima.duckdb"
BASE_AURA  <- "~/MEGA/CIE things/ducks/aurassure.duckdb"
INMET_RIO  <- c("VILA MILITAR", "RIO DE JANEIRO - JACAREPAGUA", "MARAMBAIA")

# Magnus (Alduchov & Eskridge 1996), T in C, RH in %
td_magnus <- function(t, ur) {
  g <- log(ur / 100) + 17.625 * t / (243.04 + t)
  243.04 * g / (17.625 - g)
}
# The same expression in SQL (DuckDB), for record-level computation
sql_td <- function(t, ur) {
  sprintf("(243.04 * (ln(%2$s / 100.0) + 17.625 * %1$s / (243.04 + %1$s)) / (17.625 - (ln(%2$s / 100.0) + 17.625 * %1$s / (243.04 + %1$s))))",
          t, ur)
}

buscar_orvalho_periodo <- function(data_ini, data_fim, estatisticas = c("median", "max", "min"),
                                   min_horas = 18, max_frac_saturada = 1,
                                   fontes = c("inmet", "redemet", "ar", "aura", "wunder")) {
  ini <- format(as.Date(data_ini), "%Y-%m-%d 00:00:00")
  fim <- format(as.Date(data_fim) + 1, "%Y-%m-%d 00:00:00")
  cat(sprintf("[orvalho] %s a %s\n", ini, fim))

  # Hourly Td -> daily statistics, shared by all sources
  diario <- function(con, horario_sql) {
    dbGetQuery(con, glue("
      WITH horario AS ({horario_sql})
      SELECT estacao_id, date_trunc('day', hora) AS data,
             COUNT(*) AS n_horas,
             MEDIAN(td) AS median_orvalho, MAX(td) AS max_orvalho, MIN(td) AS min_orvalho,
             AVG(saturada) AS frac_saturada
      FROM horario GROUP BY estacao_id, date_trunc('day', hora)"))
  }
  conta <- function(con, sql) dbGetQuery(con, sql)[[1]]
  relata <- function(nome, total, validos) {
    cat(sprintf("[orvalho] %-8s %9d registros, %9d validos (%.1f%% removidos pelos filtros)\n",
                nome, total, validos, 100 * (1 - validos / max(total, 1))))
  }
  qc_td <- "td BETWEEN 0 AND 30 AND (t IS NULL OR td <= t + 0.5)"

  saida <- list()
  con <- dbConnect(duckdb(), dbdir = path.expand(BASE_CLIMA), read_only = TRUE)
  coords <- dbGetQuery(con, "SELECT * FROM tb_estacoes_all")

  if ("inmet" %in% fontes) {
    ids <- coords$estacao_id[coords$estacao_nome %in% INMET_RIO]
    lista <- paste0("'", ids, "'", collapse = ", ")
    base <- glue("SELECT estacao_id, dt_registro AS dt, temperatura_bulbo AS t, umidade AS ur,
                         temperatura_orvalho AS td
                  FROM tb_clima_inmet
                  WHERE dt_registro >= '{ini}' AND dt_registro < '{fim}' AND estacao_id IN ({lista})")
    total <- conta(con, glue("SELECT COUNT(*) FROM ({base})"))
    filtrado <- glue("SELECT * FROM ({base}) WHERE {qc_td} AND (ur IS NULL OR ur BETWEEN 10 AND 100)")
    relata("INMET", total, conta(con, glue("SELECT COUNT(*) FROM ({filtrado})")))
    saida$inmet <- diario(con, glue("
      SELECT estacao_id, date_trunc('hour', dt) AS hora, MEDIAN(td) AS td,
             CAST(NULL AS DOUBLE) AS saturada
      FROM ({filtrado}) GROUP BY estacao_id, date_trunc('hour', dt)")) |>
      mutate(origem = "inmet", td_fonte = "medido")
  }

  if ("redemet" %in% fontes) {
    base <- glue("SELECT cod_aeroporto AS estacao_id, dt_referencia AS dt, temperatura AS t,
                         umidade AS ur, ponto_orvalho AS td
                  FROM tb_clima_redemet
                  WHERE dt_referencia >= '{ini}' AND dt_referencia < '{fim}'")
    total <- conta(con, glue("SELECT COUNT(*) FROM ({base})"))
    filtrado <- glue("SELECT * FROM ({base}) WHERE {qc_td} AND (ur IS NULL OR ur BETWEEN 10 AND 100)")
    relata("REDEMET", total, conta(con, glue("SELECT COUNT(*) FROM ({filtrado})")))
    saida$redemet <- diario(con, glue("
      SELECT estacao_id, date_trunc('hour', dt) AS hora, MEDIAN(td) AS td,
             CAST(NULL AS DOUBLE) AS saturada
      FROM ({filtrado}) GROUP BY estacao_id, date_trunc('hour', dt)")) |>
      mutate(origem = "redemet", td_fonte = "medido")
  }

  if ("ar" %in% fontes) {
    base <- glue("SELECT CAST(estacao_id AS VARCHAR) AS estacao_id, dt_registro AS dt,
                         temperatura AS t, umidade AS ur
                  FROM tb_clima_ar
                  WHERE dt_registro >= '{ini}' AND dt_registro < '{fim}'")
    total <- conta(con, glue("SELECT COUNT(*) FROM ({base})"))
    calc <- glue("SELECT *, {sql_td('t', 'ur')} AS td FROM ({base})
                  WHERE ur BETWEEN 10 AND 100 AND t BETWEEN 5 AND 45")
    filtrado <- glue("SELECT * FROM ({calc}) WHERE {qc_td}")
    relata("AlertaRio", total, conta(con, glue("SELECT COUNT(*) FROM ({filtrado})")))
    saida$ar <- diario(con, glue("
      SELECT estacao_id, date_trunc('hour', dt) AS hora, MEDIAN(td) AS td,
             CASE WHEN MEDIAN(ur) >= 99.5 THEN 1.0 ELSE 0.0 END AS saturada
      FROM ({filtrado}) GROUP BY estacao_id, date_trunc('hour', dt)")) |>
      mutate(origem = "ar", td_fonte = "calculado")
  }
  dbDisconnect(con, shutdown = TRUE)

  coords_of <- coords |>
    transmute(estacao_id = as.character(estacao_id), latitude, longitude) |>
    distinct(estacao_id, .keep_all = TRUE)

  if ("aura" %in% fontes) {
    con2 <- dbConnect(duckdb(), dbdir = path.expand(BASE_AURA), read_only = TRUE)
    base <- glue("SELECT CAST(CAST(thing_id AS BIGINT) AS VARCHAR) AS estacao_id, time AS dt,
                         temp AS t, humid AS ur
                  FROM tb_registers WHERE time >= '{ini}' AND time < '{fim}'")
    total <- conta(con2, glue("SELECT COUNT(*) FROM ({base})"))
    calc <- glue("SELECT *, {sql_td('t', 'ur')} AS td FROM ({base})
                  WHERE ur BETWEEN 10 AND 100 AND t BETWEEN 5 AND 45")
    filtrado <- glue("SELECT * FROM ({calc}) WHERE {qc_td}")
    relata("Aurassure", total, conta(con2, glue("SELECT COUNT(*) FROM ({filtrado})")))
    if (total > 0) {
      saida$aura <- diario(con2, glue("
        SELECT estacao_id, date_trunc('hour', dt) AS hora, MEDIAN(td) AS td,
               CASE WHEN MEDIAN(ur) >= 99.5 THEN 1.0 ELSE 0.0 END AS saturada
        FROM ({filtrado}) GROUP BY estacao_id, date_trunc('hour', dt)")) |>
        mutate(origem = "aura", td_fonte = "calculado")
      # tb_things has latitude and longitude swapped (as in the temperature query)
      coords_aura <- dbGetQuery(con2, "SELECT thing_id, lat, long FROM tb_things") |>
        transmute(estacao_id = as.character(as.integer(thing_id)), latitude = long, longitude = lat)
      coords_of <- bind_rows(coords_of, coords_aura)
    }
    dbDisconnect(con2, shutdown = TRUE)
  }

  if ("wunder" %in% fontes) {
    wu <- buscar_wu_orvalho(as.Date(data_ini), as.Date(data_fim))
    if (nrow(wu)) {
      saida$wunder <- wu |> mutate(origem = "wunder", td_fonte = "calculado")
      coords_of <- bind_rows(coords_of, wu |> distinct(estacao_id, latitude, longitude))
      saida$wunder <- saida$wunder |> select(-latitude, -longitude)
    }
  }

  td <- bind_rows(saida) |>
    mutate(estacao_id = as.character(estacao_id), data = as.Date(data)) |>
    left_join(coords_of, by = "estacao_id")

  n0 <- nrow(td)
  td <- td |> filter(n_horas >= min_horas, is.na(frac_saturada) | frac_saturada <= max_frac_saturada)
  cat(sprintf("[orvalho] %d estacao-dias; %d removidos por ter menos de %d horas%s\n",
              nrow(td), n0 - nrow(td), min_horas,
              if (max_frac_saturada < 1) sprintf(" ou mais de %.0f%% de horas saturadas", 100 * max_frac_saturada) else ""))
  sem_coord <- unique(td$estacao_id[is.na(td$latitude)])
  if (length(sem_coord)) cat("[orvalho] sem coordenadas (ficam de fora no bruce):", paste(sem_coord, collapse = ", "), "\n")

  manter <- c("estacao_id", "data", "origem", "td_fonte", "latitude", "longitude", "n_horas",
              "frac_saturada", paste0(estatisticas, "_orvalho"), "median_orvalho_wu")
  td[, intersect(manter, names(td))]
}

# ---- Wunderground ---------------------------------------------------------
STATIONS_WU <- c("IRIODE81","IRIODE113","IRIODE165","IRIODE112","IRIODE109","IRIODE159",
                 "IRIODE96","IRIODE83","IRIODE84","IRIODEJA96","IRIODE130","IRIODE157",
                 "IRIODE89","IRIODE150","IRIODE135","IRIODE114","IRIODE122","IRIODE149",
                 "IRIODE144","IRIODE85","IRIODE108")

# One station-day from the API: every observation with T, RH and WU's own Td
wu_observacoes <- function(station_id, dia, chave) {
  req <- request("https://api.weather.com/v2/pws/history/all") |>
    req_url_query(stationId = station_id, format = "json", units = "m",
                  date = format(dia, "%Y%m%d"), apiKey = chave) |>
    req_headers("Accept-Encoding" = "gzip") |>
    req_timeout(15)
  resp <- tryCatch(req_perform(req), httr2_http = function(e) e$resp, error = function(e) NULL)
  if (is.null(resp) || resp_status(resp) == 204 || resp_status(resp) >= 400) return(NULL)
  obs <- resp_body_json(resp)$observations
  if (is.null(obs) || !length(obs)) return(NULL)
  num <- function(x) if (is.null(x)) NA_real_ else as.numeric(x)
  data.frame(
    obs_time = as.POSIXct(vapply(obs, function(o) o$obsTimeUtc %||% NA_character_, ""),
                          format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    t = vapply(obs, function(o) num(o$metric$tempAvg), 0),
    ur = vapply(obs, function(o) num(o$humidityAvg), 0),
    td_wu = vapply(obs, function(o) num(o$metric$dewptAvg), 0),
    lat = vapply(obs, function(o) num(o$lat), 0),
    lon = vapply(obs, function(o) num(o$lon), 0))
}

# Same filters and two stages as the other sources; one row per station-day
wu_resumo_dia <- function(o, station_id, dia) {
  o <- o[!is.na(o$t) & !is.na(o$ur) & o$ur >= 10 & o$ur <= 100 & o$t >= 5 & o$t <= 45, ]
  if (!nrow(o)) return(NULL)
  o$td <- td_magnus(o$t, o$ur)
  o <- o[o$td >= 0 & o$td <= 30 & o$td <= o$t + 0.5, ]
  if (!nrow(o)) return(NULL)
  o$hora <- format(o$obs_time, "%Y-%m-%d %H")
  h <- do.call(rbind, lapply(split(o, o$hora), function(x) data.frame(
    td = stats::median(x$td), td_wu = stats::median(x$td_wu, na.rm = TRUE),
    saturada = as.numeric(stats::median(x$ur) >= 99.5))))
  data.frame(estacao_id = station_id, data = as.Date(dia), n_horas = nrow(h),
             median_orvalho = stats::median(h$td), max_orvalho = max(h$td), min_orvalho = min(h$td),
             frac_saturada = mean(h$saturada), median_orvalho_wu = stats::median(h$td_wu, na.rm = TRUE),
             latitude = o$lat[1], longitude = o$lon[1])
}

buscar_wu_orvalho <- function(ini, fim, estacoes = STATIONS_WU,
                              cache = file.path(dirname(path.expand(
                                if (exists("WU_CACHE")) WU_CACHE else "wu_cache_v2.duckdb")),
                                "wu_orvalho_cache.duckdb")) {
  con <- dbConnect(duckdb(), dbdir = cache)
  on.exit(dbDisconnect(con, shutdown = TRUE))
  dbExecute(con, "CREATE TABLE IF NOT EXISTS wu_orvalho (estacao_id VARCHAR, data DATE,
    n_horas INTEGER, median_orvalho DOUBLE, max_orvalho DOUBLE, min_orvalho DOUBLE,
    frac_saturada DOUBLE, median_orvalho_wu DOUBLE, latitude DOUBLE, longitude DOUBLE)")
  # station-days the API returned nothing for, so they are not asked again
  dbExecute(con, "CREATE TABLE IF NOT EXISTS wu_orvalho_vazio (estacao_id VARCHAR, data DATE)")
  dias <- seq(ini, fim, by = "day")
  todos <- expand.grid(estacao_id = estacoes, data = dias, stringsAsFactors = FALSE)
  feitos <- rbind(dbGetQuery(con, "SELECT estacao_id, data FROM wu_orvalho"),
                  dbGetQuery(con, "SELECT estacao_id, data FROM wu_orvalho_vazio"))
  chave_feito <- paste(feitos$estacao_id, as.Date(feitos$data))
  faltam <- todos[!paste(todos$estacao_id, todos$data) %in% chave_feito, ]
  chave <- Sys.getenv("WU_API_KEY")
  cat(sprintf("[orvalho] Wunderground: %d estacao-dias em cache, %d a buscar%s\n",
              nrow(todos) - nrow(faltam), nrow(faltam),
              if (nrow(faltam) && !nzchar(chave)) " (sem WU_API_KEY: so o cache)" else ""))
  if (nrow(faltam) && nzchar(chave)) {
    for (i in seq_len(nrow(faltam))) {
      if (i %% 50 == 0 || i == nrow(faltam)) cat(sprintf("  [orvalho wu] %d/%d\n", i, nrow(faltam)))
      Sys.sleep(0.2)
      o <- wu_observacoes(faltam$estacao_id[i], faltam$data[i], chave)
      r <- if (is.null(o)) NULL else wu_resumo_dia(o, faltam$estacao_id[i], faltam$data[i])
      if (is.null(r)) {
        dbWriteTable(con, "wu_orvalho_vazio", data.frame(estacao_id = faltam$estacao_id[i],
                                                         data = as.Date(faltam$data[i])), append = TRUE)
      } else {
        dbWriteTable(con, "wu_orvalho", r, append = TRUE)
      }
    }
  }
  out <- dbGetQuery(con, glue("SELECT * FROM wu_orvalho WHERE data >= '{ini}' AND data <= '{fim}'"))
  out$data <- as.Date(out$data)
  out
}
