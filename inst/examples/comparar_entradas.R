# ==================================================================
# comparar_entradas.R
# Compares the INPUTS of the article's fit (saved .rds) with those of the
# bruce fit from rio_adapter.R, without refitting anything: which stations
# and days are in each, whether the observed values agree, whether the
# covariates at stations and at ERA5 cells agree, and the mesh size.
# Run from ~/scratch/spde_paper. Takes seconds.
# ==================================================================

library(sf)
library(dplyr)

VAR <- "median"; MES <- "2025-01"
orig <- readRDS(sprintf("resultados_9meses_median_max/pacote_mensal_independente_%s_%s.rds", VAR, MES))
novo <- readRDS(sprintf("resultados_bruce/bruce_%s_%s.rds", VAR, MES))
d    <- novo$fit$data

# ---- 1. Stations and observations ------------------------------------------
est <- st_drop_geometry(orig$est_all)
est$source <- unname(c(ofic = "official", aura = "aurassure", wunder = "wunder")[est$origem])
est$source[est$estacao_id %in% orig$ESTACOES_AEROPORTO] <- "airport"
a <- data.frame(station_id = as.character(est$estacao_id), date = as.Date(est$data),
                source = est$source, valor_artigo = est[[orig$col_estacao]])
b <- data.frame(station_id = d$obs$station_id, date = d$obs$date,
                source = d$obs$source, valor_bruce = d$obs$value)

cat("=== Estacoes e observacoes por fonte ===\n")
resumo <- full_join(
  a %>% group_by(source) %>% summarise(est_artigo = n_distinct(station_id), obs_artigo = n()),
  b %>% group_by(source) %>% summarise(est_bruce = n_distinct(station_id), obs_bruce = n()),
  by = "source")
print(as.data.frame(resumo), row.names = FALSE)

so_artigo <- setdiff(unique(a$station_id), unique(b$station_id))
so_bruce  <- setdiff(unique(b$station_id), unique(a$station_id))
cat(sprintf("\nEstacoes so no artigo (%d): %s\n", length(so_artigo), paste(head(so_artigo, 15), collapse = ", ")))
cat(sprintf("Estacoes so no bruce  (%d): %s\n", length(so_bruce), paste(head(so_bruce, 15), collapse = ", ")))

ab <- inner_join(a, b, by = c("station_id", "date"), suffix = c("_artigo", "_bruce"))
ab$dif <- ab$valor_bruce - ab$valor_artigo
cat(sprintf("\nObservacoes em comum: %d (artigo %d, bruce %d)\n", nrow(ab), nrow(a), nrow(b)))
cat("Diferenca de valor nas observacoes em comum, por fonte:\n")
print(as.data.frame(ab %>% group_by(source = source_artigo) %>%
  summarise(n = n(), n_diferentes = sum(abs(dif) > 1e-6), max_abs = max(abs(dif)),
            media = mean(dif))), row.names = FALSE, digits = 3)
mudou_fonte <- ab %>% filter(source_artigo != source_bruce) %>% distinct(station_id, source_artigo, source_bruce)
if (nrow(mudou_fonte)) { cat("\nEstacoes classificadas em fontes diferentes:\n"); print(mudou_fonte) }

# ---- 2. Covariates at stations (raw values) --------------------------------
sc  <- d$covariates$scaling
raw <- function(z, nm) z * sc$sd[sc$name == nm] + sc$mean[sc$name == nm]
cov_b <- d$covariates$station %>%
  transmute(station_id, elev_bruce = raw(elev, "elev"), urb_bruce = raw(urb, "urb"), svf_bruce = raw(svf, "svf"))
cov_a <- est %>% transmute(station_id = as.character(estacao_id), elev_artigo = elev,
                           urb_artigo = frac_urb, svf_artigo = svf) %>% distinct(station_id, .keep_all = TRUE)
cv <- inner_join(cov_a, cov_b, by = "station_id")
cat("\n=== Covariaveis nas estacoes (valores brutos) ===\n")
for (v in c("elev", "urb", "svf")) {
  x <- cv[[paste0(v, "_artigo")]]; y <- cv[[paste0(v, "_bruce")]]; ok <- !is.na(x) & !is.na(y)
  cat(sprintf("%-4s  n=%d  cor=%.4f  max|dif|=%.3f  estacoes com |dif|>1%% da amplitude: %d\n",
              v, sum(ok), cor(x[ok], y[ok]), max(abs(x[ok] - y[ok])),
              sum(abs(x[ok] - y[ok]) > 0.01 * diff(range(x[ok])))))
}

# ---- 3. Background (ERA5) ---------------------------------------------------
bg_a <- st_drop_geometry(orig$d_all)
cat(sprintf("\n=== ERA5 ===\nCelulas x dias: artigo %d | bruce %d\n", nrow(bg_a), nrow(d$background)))
cat(sprintf("Celulas distintas: artigo %d | bruce %d\n",
            nrow(distinct(bg_a, x, y)), nrow(distinct(d$background, x, y))))
cb <- d$covariates$background %>% transmute(x, y, elev_b = raw(elev, "elev"), urb_b = raw(urb, "urb"))
ca <- bg_a %>% distinct(x, y, .keep_all = TRUE) %>% transmute(x, y, elev_a = elev, urb_a = frac_urb)
m  <- inner_join(ca, cb, by = c("x", "y"))
cat(sprintf("Celulas casadas por coordenada: %d\n", nrow(m)))
if (nrow(m)) {
  cat(sprintf("elev  max|dif| = %.1f m | urb max|dif| = %.3f | celulas com urb diferente: %d\n",
              max(abs(m$elev_a - m$elev_b)), max(abs(m$urb_a - m$urb_b)),
              sum(abs(m$urb_a - m$urb_b) > 0.01)))
}

# ---- 4. Mesh ------------------------------------------------------------------
cat(sprintf("\nVertices da malha: artigo %d | bruce %d\n", orig$mesh$n, novo$fit$mesh$n))
