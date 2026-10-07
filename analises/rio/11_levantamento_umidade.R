# ==================================================================
# 11_levantamento_umidade.R
# Which networks measure humidity, and how good is it? First step before
# modelling dew point with bruce. Reads the DuckDB databases (read only),
# fits nothing.
#
# For each table used by query_meteorj_periodo.R, finds columns that look
# like humidity or dew point and, for one month, reports: stations with
# data, records, fraction missing, the distribution of values, and signs of
# trouble (relative humidity above 100%, stuck at 99-100%, dew point above
# temperature). Wunderground is not here: its cache keeps temperature only.
#
# Rode de dentro de analises/rio:  Rscript 11_levantamento_umidade.R
# ==================================================================

library(DBI)
library(duckdb)

INI <- "2025-01-01"; FIM <- "2025-02-01"
BASES <- list(
  clima = list(arquivo = "~/MEGA/CIE things/ducks/clima.duckdb",
               tabelas = c(tb_clima_inmet = "dt_registro", tb_clima_ar = "dt_registro",
                           tb_clima_redemet = "dt_referencia"),
               estacao = c(tb_clima_inmet = "estacao_id", tb_clima_ar = "estacao_id",
                           tb_clima_redemet = "cod_aeroporto")),
  aurassure = list(arquivo = "~/MEGA/CIE things/ducks/aurassure.duckdb",
                   tabelas = c(tb_registers = "time"), estacao = c(tb_registers = "thing_id"))
)
padrao <- "umid|humid|orvalho|dew|(^|_)(ur|rh)($|_)|ponto"
pad_temp <- "^temp|temperatura"

linhas <- list()
for (b in names(BASES)) {
  cfg <- BASES[[b]]
  con <- dbConnect(duckdb(), dbdir = path.expand(cfg$arquivo), read_only = TRUE)
  for (tab in names(cfg$tabelas)) {
    cols <- dbListFields(con, tab)
    cat(sprintf("\n== %s.%s ==\nColunas: %s\n", b, tab, paste(cols, collapse = ", ")))
    hum <- grep(padrao, cols, value = TRUE, ignore.case = TRUE)
    if (!length(hum)) { cat("  nenhuma coluna de umidade ou ponto de orvalho\n"); next }
    tcol <- cfg$tabelas[[tab]]; est <- cfg$estacao[[tab]]
    temp <- grep(pad_temp, cols, value = TRUE, ignore.case = TRUE)[1]
    for (h in hum) {
      q <- sprintf(paste0(
        "SELECT COUNT(DISTINCT %s) AS estacoes, COUNT(*) AS registros, ",
        "AVG(CASE WHEN %s IS NULL THEN 1.0 ELSE 0 END) AS frac_vazio, ",
        "MIN(%s) AS min, QUANTILE_CONT(%s, 0.01) AS p01, MEDIAN(%s) AS mediana, ",
        "QUANTILE_CONT(%s, 0.99) AS p99, MAX(%s) AS max, ",
        "AVG(CASE WHEN %s > 100 THEN 1.0 ELSE 0 END) AS frac_acima_100, ",
        "AVG(CASE WHEN %s >= 99 AND %s <= 100 THEN 1.0 ELSE 0 END) AS frac_99_100%s ",
        "FROM %s WHERE %s >= '%s' AND %s < '%s'"),
        est, h, h, h, h, h, h, h, h, h,
        if (!is.na(temp)) sprintf(", AVG(CASE WHEN %s > %s THEN 1.0 ELSE 0 END) AS frac_acima_da_temp", h, temp) else "",
        tab, tcol, INI, tcol, FIM)
      r <- tryCatch(dbGetQuery(con, q), error = function(e) { message("  ", h, ": ", conditionMessage(e)); NULL })
      if (is.null(r)) next
      linhas[[paste(b, tab, h)]] <- data.frame(base = b, tabela = tab, coluna = h, r)
    }
  }
  dbDisconnect(con, shutdown = TRUE)
}

if (length(linhas)) {
  tab <- do.call(rbind, linhas); rownames(tab) <- NULL
  num <- vapply(tab, is.numeric, TRUE)
  tab[num] <- lapply(tab[num], function(x) round(x, 3))
  cat(sprintf("\n=== Colunas de umidade / ponto de orvalho, %s a %s ===\n", INI, FIM))
  print(tab, row.names = FALSE)
  cat("frac_acima_100 e frac_99_100: so fazem sentido para umidade relativa (%).\n")
  cat("frac_acima_da_temp: so faz sentido para ponto de orvalho (deve ser ~0).\n")
  write.csv(tab, "levantamento_umidade.csv", row.names = FALSE)
} else {
  cat("\nNenhuma coluna de umidade encontrada nas tabelas.\n")
}
