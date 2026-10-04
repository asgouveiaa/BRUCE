# ==================================================================
# query_meteorj_periodo_v2.R
# Versao de TESTE: corrige o viés de cadencia de amostragem
# identificado (tb_clima_ar ~96 leituras/dia, Wunderground ~288/dia,
# vs. INMET/REDEMET ~24/dia) -- reamostra para HORARIO primeiro
# (MEDIAN/MAX por hora), so DEPOIS calcula mediana/maximo do dia a
# partir dessas ~24 estatisticas horarias, igualando o criterio
# entre todas as fontes.
#
# NAO mexe no TODO 2 (REDEMET usa sensacao, nao temperatura) --
# isolado deliberadamente, para nao misturar duas correcoes na mesma
# comparacao e perder a capacidade de atribuir causa.
#
# INMET fica inalterado -- MAX(temperatura_max_hora) ja e derivado de
# 24 maximos horarios nativos, nao sofre desse problema.
#
# v3 (out/2026): temperatura MINIMA opcional.
#   buscar_meteo_periodo_bulk(ini, fim, estatisticas = c("median", "max", "min"))
#   * Default c("median", "max"): saida IDENTICA a v2 (mesmas linhas, mesmas
#     colunas, mesmos valores) -- reprodutibilidade dos resultados do artigo.
#   * Com "min": acrescenta min_temperatura, calculada como a mediana/maximo
#     (horario primeiro, dia depois): MIN dos minimos horarios. O filtro
#     fisico da minima so anula a propria minima (vira NA); nao remove a
#     linha, para nao alterar mediana e maxima.
#   * Wunderground: o cache ganha a coluna min_temperatura. Dias ja em cache
#     sem minima sao buscados de novo na API, mas so a minima e gravada --
#     mediana e maxima ja cacheadas nao mudam.
# ==================================================================

library(tidyverse)
library(duckdb)
library(DBI)
library(glue)
library(lubridate)
library(httr2)

# ---- Wunderground: agora preserva o timestamp da observacao, para
# poder reamostrar por hora antes de agregar o dia -------------------------
buscar_wu_dia <- function(station_id, dia) {
  if (!nzchar(api_wunder)) return(NULL)   # sem chave: so o que ja estiver no cache
  req <- request("https://api.weather.com/v2/pws/history/all") %>%
    req_url_query(stationId = station_id, format = "json", units = "m",
                  date = format(dia, "%Y%m%d"), apiKey = api_wunder) %>%
    req_headers("Accept-Encoding" = "gzip") %>%
    req_timeout(15)   # evita travamento silencioso se a rede travar nesta chamada
  
  resp <- tryCatch(req_perform(req), httr2_http = function(e) e$resp, error = function(e) NULL)
  if (is.null(resp) || resp_status(resp) == 204 || resp_status(resp) >= 400) return(NULL)
  
  obs <- resp_body_json(resp)$observations
  if (is.null(obs) || length(obs) == 0) return(NULL)
  
  map_dfr(obs, function(o) tibble(
    estacao_id = o$stationID %||% station_id,
    obs_time   = as_datetime(o$obsTimeUtc %||% NA, tz = "UTC"),
    temp = o$metric$tempAvg %||% NA_real_,
    lat = o$lat %||% NA_real_, lon = o$lon %||% NA_real_
  ))
}

buscar_wu_periodo <- function(stations_wu, dias,
                              cache_path = if (exists("WU_CACHE")) WU_CACHE else "wu_cache_v2.duckdb",
                              incluir_min = FALSE) {
  
  # Conexao com retry -- quando rodando varias janelas em paralelo
  # (rodar_piloto_sazonal.sh), mais de um processo pode tentar abrir o
  # mesmo wu_cache_v2.duckdb ao mesmo tempo. DuckDB nao permite duas
  # conexoes de escrita simultaneas no mesmo arquivo (lock exclusivo).
  # Em vez de falhar na primeira tentativa, espera um pouco e tenta de
  # novo -- o lock e sempre breve (so durante leitura/escrita do cache).
  #
  # NOTA: essa blindagem ja tinha sido adicionada antes, mas foi perdida
  # quando o arquivo foi reescrito para a versao com reamostragem
  # horaria (v2) -- reintroduzida aqui apos o mesmo erro reaparecer.
  conectar_cache_com_retry <- function(cache_path, max_tentativas = 10) {
    for (tentativa in seq_len(max_tentativas)) {
      con <- tryCatch(
        dbConnect(duckdb(), dbdir = cache_path),
        error = function(e) {
          if (grepl("lock", conditionMessage(e), ignore.case = TRUE)) {
            espera <- runif(1, 2, 8) * tentativa   # backoff crescente
            cat(sprintf("[wu_cache] lock ocupado, tentativa %d/%d, esperando %.1fs...\n",
                        tentativa, max_tentativas, espera))
            Sys.sleep(espera)
            return(NULL)
          } else {
            stop(e)   # erro diferente de lock -- nao adianta tentar de novo
          }
        }
      )
      if (!is.null(con)) return(con)
    }
    stop(sprintf("[wu_cache] nao foi possivel conectar a %s apos %d tentativas (lock persistente)",
                 cache_path, max_tentativas))
  }
  
  con_cache <- conectar_cache_com_retry(cache_path)
  dbExecute(con_cache, "CREATE TABLE IF NOT EXISTS wu_cache (
    estacao_id VARCHAR, dia DATE, median_temperatura DOUBLE,
    max_temperatura DOUBLE, latitude DOUBLE, longitude DOUBLE)")
  # v3: coluna da minima (caches antigos ganham a coluna, com NULL)
  dbExecute(con_cache, "ALTER TABLE wu_cache ADD COLUMN IF NOT EXISTS min_temperatura DOUBLE")
  
  vazio <- tibble(estacao_id = character(), dia = as_date(character()))
  cache_todas <- tryCatch(dbGetQuery(con_cache, "SELECT estacao_id, dia FROM wu_cache"),
                          error = function(e) vazio)
  # Com incluir_min, um dia so conta como "em cache" se ja tiver a minima
  cache_existente <- if (incluir_min) {
    tryCatch(dbGetQuery(con_cache, "SELECT estacao_id, dia FROM wu_cache WHERE min_temperatura IS NOT NULL"),
             error = function(e) vazio)
  } else cache_todas
  
  combinacoes <- expand.grid(estacao_id = stations_wu, dia = dias, stringsAsFactors = FALSE) %>%
    mutate(dia = as_date(dia))
  faltantes <- combinacoes %>% anti_join(cache_existente, by = c("estacao_id", "dia"))
  
  cat(sprintf("[wu v2] %d combinacoes ja em cache, buscando %d novas via API\n",
              nrow(combinacoes) - nrow(faltantes), nrow(faltantes)))
  
  if (nrow(faltantes) > 0) {
    n_total_faltantes <- nrow(faltantes)
    novos <- map_dfr(seq_len(n_total_faltantes), function(i) {
      estacao_id <- faltantes$estacao_id[i]
      dia        <- faltantes$dia[i]
      
      if (i %% 50 == 0 || i == n_total_faltantes) {
        cat(sprintf("  [wu v2] progresso: %d/%d combinacoes buscadas\n", i, n_total_faltantes))
      }
      
      Sys.sleep(0.2)
      obs_dia <- buscar_wu_dia(estacao_id, dia)
      if (is.null(obs_dia) || nrow(obs_dia) == 0) return(NULL)
      
      # REAMOSTRAGEM: agrega por HORA primeiro, so depois pelo dia --
      # iguala a cadencia efetiva ao das fontes oficiais (~24 pontos/dia)
      obs_horaria <- obs_dia %>%
        mutate(hora = floor_date(obs_time, "hour")) %>%
        group_by(hora) %>%
        summarise(med_hora = median(temp, na.rm = TRUE), max_hora = max(temp, na.rm = TRUE),
                  min_hora = min(temp, na.rm = TRUE), .groups = "drop")
      
      tibble(
        estacao_id = first(obs_dia$estacao_id), dia = dia,
        median_temperatura = median(obs_horaria$med_hora, na.rm = TRUE),
        max_temperatura    = max(obs_horaria$max_hora, na.rm = TRUE),
        latitude = first(obs_dia$lat), longitude = first(obs_dia$lon),
        min_temperatura    = min(obs_horaria$min_hora, na.rm = TRUE)
      )
    })
    
    if (nrow(novos) > 0) {
      novos <- novos %>% mutate(dia = as_date(dia))
      ja_tinha <- novos %>% semi_join(cache_todas %>% mutate(dia = as_date(dia)),
                                      by = c("estacao_id", "dia"))
      inedito  <- novos %>% anti_join(cache_todas %>% mutate(dia = as_date(dia)),
                                      by = c("estacao_id", "dia"))
      if (nrow(inedito) > 0) {
        dbWriteTable(con_cache, "wu_cache",
                     inedito %>% select(estacao_id, dia, median_temperatura, max_temperatura,
                                        latitude, longitude, min_temperatura),
                     append = TRUE)
      }
      if (nrow(ja_tinha) > 0) {
        # so preenche a minima -- mediana/maxima ja cacheadas ficam como estavam
        duckdb::duckdb_register(con_cache, "wu_min_tmp",
                                ja_tinha %>% select(estacao_id, dia, min_temperatura))
        dbExecute(con_cache, "
          UPDATE wu_cache SET min_temperatura = t.min_temperatura
          FROM wu_min_tmp t
          WHERE wu_cache.estacao_id = t.estacao_id AND wu_cache.dia = t.dia")
        duckdb::duckdb_unregister(con_cache, "wu_min_tmp")
      }
    }
  }
  
  resultado <- dbGetQuery(con_cache, glue_sql("
    SELECT * FROM wu_cache WHERE dia IN ({dias*}) AND estacao_id IN ({stations_wu*})
  ", .con = con_cache))
  
  dbDisconnect(con_cache, shutdown = TRUE)
  if (!incluir_min) resultado <- resultado %>% select(-any_of("min_temperatura"))
  resultado
}

stations_wu <- c("IRIODE81","IRIODE113","IRIODE165","IRIODE112","IRIODE109","IRIODE159",
                 "IRIODE96","IRIODE83","IRIODE84","IRIODEJA96","IRIODE130","IRIODE157",
                 "IRIODE89","IRIODE150","IRIODE135","IRIODE114","IRIODE122","IRIODE149",
                 "IRIODE144","IRIODE85","IRIODE108")
api_wunder <- Sys.getenv("WU_API_KEY")   # chave em ~/.Renviron (WU_API_KEY=...), fora do codigo

# ---- Fontes oficiais + Aurassure: tb_clima_ar corrigido (2 estagios) ---------
buscar_meteo_periodo_bulk <- function(data_ini, data_fim, estatisticas = c("median", "max")) {
  
  if (!all(estatisticas %in% c("median", "max", "min"))) {
    stop("`estatisticas` aceita apenas 'median', 'max' e 'min'")
  }
  incluir_min <- "min" %in% estatisticas
  
  ini_dt <- format(as_date(data_ini), "%Y-%m-%d 00:00:00")
  fim_dt <- format(as_date(data_fim) + days(1), "%Y-%m-%d 00:00:00")
  
  cat(sprintf("[meteo_rj_periodo_v2] Puxando dados de %s ate %s\n", ini_dt, fim_dt))
  
  con <- dbConnect(duckdb(), dbdir = "~/MEGA/CIE things/ducks/clima.duckdb", read_only = TRUE)
  
  # Minima do INMET: usa o minimo horario nativo (espelho do MAX(temperatura_max_hora))
  # se a coluna existir; se nao, cai no MIN da temperatura de bulbo.
  campos_inmet <- dbListFields(con, "tb_clima_inmet")
  expr_min_inmet <- if ("temperatura_min_hora" %in% campos_inmet) {
    DBI::SQL("MIN(temperatura_min_hora)")
  } else {
    if (incluir_min) cat("[meteo_rj_periodo_v3] INMET sem temperatura_min_hora -- usando MIN(temperatura_bulbo)\n")
    DBI::SQL("MIN(temperatura_bulbo)")
  }
  
  meteo_inmet <- dbGetQuery(con, glue_sql("
    SELECT estacao_id, date_trunc('day', dt_registro) AS data,
           MEDIAN(temperatura_bulbo)       AS median_temperatura,
           MAX(temperatura_max_hora) AS max_temperatura,
           {expr_min_inmet}          AS min_temperatura
    FROM tb_clima_inmet
    WHERE dt_registro >= {ini_dt} AND dt_registro < {fim_dt}
    GROUP BY estacao_id, date_trunc('day', dt_registro)
  ", .con = con))
  
  # CORRIGIDO: 2 estagios -- horario primeiro, dia depois
  meteo_ar <- dbGetQuery(con, glue_sql("
    WITH horario AS (
      SELECT estacao_id, date_trunc('hour', dt_registro) AS hora,
             MEDIAN(temperatura) AS med_hora, MAX(temperatura) AS max_hora,
             MIN(temperatura) AS min_hora
      FROM tb_clima_ar
      WHERE dt_registro >= {ini_dt} AND dt_registro < {fim_dt}
      GROUP BY estacao_id, date_trunc('hour', dt_registro)
    )
    SELECT estacao_id, date_trunc('day', hora) AS data,
           MEDIAN(med_hora) AS median_temperatura,
           MAX(max_hora)    AS max_temperatura,
           MIN(min_hora)    AS min_temperatura
    FROM horario
    GROUP BY estacao_id, date_trunc('day', hora)
  ", .con = con))
  
  # NAO alterado (ver nota no cabecalho -- TODO 2 isolado desta comparacao)
  meteo_redemet <- dbGetQuery(con, glue_sql("
    SELECT cod_aeroporto, date_trunc('day', dt_referencia) AS data,
           MEDIAN(temperatura) AS median_temperatura,
           MAX(temperatura)    AS max_temperatura,
           MIN(temperatura)    AS min_temperatura
    FROM tb_clima_redemet
    WHERE dt_referencia >= {ini_dt} AND dt_referencia < {fim_dt}
    GROUP BY cod_aeroporto, date_trunc('day', dt_referencia)
  ", .con = con))
  
  meteo_coords <- dbGetQuery(con, "SELECT * FROM tb_estacoes_all")
  dbDisconnect(con, shutdown = TRUE)
  
  meteo_redemet <- meteo_redemet %>%
    left_join(meteo_coords, by = c("cod_aeroporto" = "estacao_id")) %>%
    transmute(estacao_id = cod_aeroporto, data, median_temperatura, max_temperatura,
              min_temperatura, latitude, longitude)
  
  meteo_ar <- meteo_ar %>%
    mutate(estacao_id = as.character(estacao_id)) %>%
    left_join(meteo_coords, by = "estacao_id") %>%
    select(estacao_id, data, median_temperatura, max_temperatura, min_temperatura, latitude, longitude)
  
  meteo_inmet <- meteo_inmet %>%
    left_join(meteo_coords, by = "estacao_id") %>%
    filter(estacao_nome %in% c("VILA MILITAR", "RIO DE JANEIRO - JACAREPAGUA", "MARAMBAIA")) %>%
    select(estacao_id, data, median_temperatura, max_temperatura, min_temperatura, latitude, longitude)
  
  meteo_ofc <- bind_rows(meteo_ar, meteo_redemet, meteo_inmet) %>% mutate(origem = "ofic")
  
  # CORRIGIDO: Aurassure tambem em 2 estagios
  con2 <- dbConnect(duckdb(), dbdir = "~/MEGA/CIE things/ducks/aurassure.duckdb", read_only = TRUE)
  
  aurassure <- dbGetQuery(con2, glue_sql("
    WITH horario AS (
      SELECT thing_id, date_trunc('hour', time) AS hora,
             MEDIAN(temp) AS med_hora, MAX(temp) AS max_hora, MIN(temp) AS min_hora
      FROM tb_registers
      WHERE time >= {ini_dt} AND time < {fim_dt}
      GROUP BY thing_id, date_trunc('hour', time)
    )
    SELECT thing_id, date_trunc('day', hora) AS dt_registro,
           MEDIAN(med_hora) AS median_temperatura,
           MAX(max_hora)    AS max_temperatura,
           MIN(min_hora)    AS min_temperatura
    FROM horario
    GROUP BY thing_id, date_trunc('day', hora)
  ", .con = con2))
  
  aurassure_coords_raw <- dbGetQuery(con2, "SELECT * FROM tb_things")
  dbDisconnect(con2, shutdown = TRUE)
  
  meteo_aura <- aurassure %>%
    rename(estacao_id = thing_id, data = dt_registro) %>%
    mutate(estacao_id = as.character(estacao_id)) %>%
    left_join(
      aurassure_coords_raw %>%
        select(estacao_id = thing_id, longitude = lat, latitude = long) %>%
        mutate(estacao_id = as.character(estacao_id)),
      by = "estacao_id"
    ) %>%
    mutate(origem = "aura")
  
  dias <- seq(as_date(data_ini), as_date(data_fim), by = "day")
  meteo_wu <- buscar_wu_periodo(stations_wu, dias, incluir_min = incluir_min)
  if (nrow(meteo_wu) > 0) {
    meteo_wu <- meteo_wu %>% rename(data = dia) %>% mutate(origem = "wunder")
  } else {
    meteo_wu <- tibble(estacao_id = character(), data = as_date(character()),
                       median_temperatura = numeric(), max_temperatura = numeric(),
                       min_temperatura = numeric(),
                       latitude = numeric(), longitude = numeric(), origem = character())
  }
  
  meteo_rj <- bind_rows(meteo_ofc, meteo_aura, meteo_wu)
  
  n_antes <- nrow(meteo_rj)
  meteo_rj <- meteo_rj %>%
    filter(
      is.na(median_temperatura) | dplyr::between(median_temperatura, 7, 45),
      is.na(max_temperatura)    | dplyr::between(max_temperatura, 7, 45)
    )
  n_removido <- n_antes - nrow(meteo_rj)
  if (n_removido > 0) cat(sprintf("[meteo_rj_periodo_v2] %d linhas removidas por filtro fisico\n", n_removido))
  
  if (incluir_min) {
    # Filtro da minima: so anula a propria minima, sem remover a linha (mediana e
    # maxima ficam identicas a v2). Tambem anula minima acima da mediana do dia
    # (impossivel fisicamente) e valores nao finitos (dia sem leituras).
    ruim <- !is.finite(meteo_rj$min_temperatura) |
      !dplyr::between(meteo_rj$min_temperatura, 7, 45) |
      (!is.na(meteo_rj$median_temperatura) & meteo_rj$min_temperatura > meteo_rj$median_temperatura)
    ruim[is.na(ruim)] <- TRUE
    n_min_na <- sum(ruim & !is.na(meteo_rj$min_temperatura))
    meteo_rj$min_temperatura[ruim] <- NA_real_
    if (n_min_na > 0) cat(sprintf("[meteo_rj_periodo_v3] %d minimas anuladas por filtro fisico\n", n_min_na))
  } else {
    meteo_rj <- meteo_rj %>% select(-any_of("min_temperatura"))
  }
  
  cat(sprintf("[meteo_rj_periodo_v2] %d linhas geradas\n", nrow(meteo_rj)))
  meteo_rj
}