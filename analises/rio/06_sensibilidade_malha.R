# ==================================================================
# 06_sensibilidade_malha.R
# Does the estimated spatial range stabilise as the mesh gets finer?
#
# Same data, covariates and priors as the article's January 2025 median fit
# (taken from the bruce fit saved by 01_ajuste_mes.R); only the mesh changes:
#   artigo      the mesh saved with the article (fmesher of that time)
#   5000/3000   the article's settings, rebuilt with today's fmesher
#   3000/1500   finer
#   1500/750    finer
#   1000/500    finest (slowest; optional, see RODAR_MAIS_FINA)
# Coarse to fine, each mesh saved as soon as it finishes (resumable).
#
# Expectation: if range and maps converge as edges shrink below ~range/5,
# report the converged value; if the range keeps moving, its weak
# identification goes beyond the prior (Limitations).
#
# Rode de dentro de analises/rio (le os caminhos de config.R).
#   Rscript 06_sensibilidade_malha.R > log_sensibilidade_malha.log 2>&1
# ==================================================================

library(bruce)
library(terra)
library(sf)
suppressMessages(library(INLA))
source("config.R")
inla.setOption(fmesher.evolution.warn = FALSE)

VAR <- "median"; MES <- "2025-01"
RODAR_MAIS_FINA <- FALSE   # TRUE adds 1500/750 and 1000/500 (many hours, lots of RAM)
# Start each fit from the hyperparameters of the 01 fit (warm start: theta is
# still learned). The 3000/1500 mesh collapsed (GroupRho = 1, field sd ~1e5)
# when started from scratch; the article's LRO region 5 was fixed the same way.
PARTIDA_QUENTE <- TRUE
# Meshes saved with a collapsed fit are refitted (the old file is kept as
# malha_<nome>_colapsada.rds and left out of the summary).
REFAZER_COLAPSADAS <- TRUE
# Memory. INLA keeps one copy of the (large) factorised precision matrix per
# configuration it evaluates in parallel: num.threads "A:B" -> A copies.
# With all 20 cores in A, the 3000/1500 mesh (3,601 vertices x 31 days) came
# close to 64 GB. A = 4 (or 2 for the finer meshes) keeps it well below.
# "eb" reproduces "ccd" (08_theta_emprestado.R, 9 windows) with fewer evaluations.
INT_STRATEGY <- "eb"
NUM_THREADS  <- "4:1"    # use "2:1" for 1500/750 and finer
N_DIAS       <- NULL     # NULL = whole month; e.g. 15 roughly halves memory and time
DIR <- file.path(RES, sprintf("sensibilidade_malha_%s%s", INT_STRATEGY,
                              if (is.null(N_DIAS)) "" else sprintf("_%dd", N_DIAS)))
dir.create(DIR, showWarnings = FALSE, recursive = TRUE)

orig <- readRDS(sprintf(file.path(ARTIGO, "pacote_mensal_independente_%s_%s.rds"), VAR, MES))
novo <- readRDS(sprintf(file.path(RES, "bruce_%s_%s.rds"), VAR, MES))
d    <- novo$fit$data
spec <- bruce_spec(range_prior = 15000, intercept_prior = c(0, 0.001), on_collapse = "warn",
                   int_strategy = INT_STRATEGY, num_threads = NUM_THREADS)
janela <- range(d$obs$date)
if (!is.null(N_DIAS)) janela[2] <- janela[1] + N_DIAS - 1

malhas <- list(
  artigo      = function() orig$mesh,
  `5000/3000` = function() bruce_mesh(d, max_edge = c(5000, 40000), cutoff = 3000, offset = c(1000, 50000)),
  `3000/1500` = function() bruce_mesh(d, max_edge = c(3000, 40000), cutoff = 1500, offset = c(1000, 50000))
)
if (RODAR_MAIS_FINA) {
  malhas$`1500/750` <- function() bruce_mesh(d, max_edge = c(1500, 40000), cutoff = 750, offset = c(1000, 50000))
  malhas$`1000/500` <- function() bruce_mesh(d, max_edge = c(1000, 40000), cutoff = 500, offset = c(1000, 50000))
}

cat("Vertices de cada malha:\n")
for (nm in names(malhas)) cat(sprintf("  %-10s %6d\n", nm, malhas[[nm]]()$n))

# 500 m grid of the article, for the map comparison
mrj <- geobr::read_municipality(code_muni = 3304557, year = 2020) |> st_transform(31983)
grid500 <- rast(vect(mrj), resolution = 500, crs = "EPSG:31983"); values(grid500) <- 1
grid500 <- mask(grid500, vect(mrj))
rs <- c(elev = "average", urb = "near", svf = "near")

pred_art <- read.csv(sprintf(file.path(ARTIGO, "predicao_mensal_independente_%s_%s.csv"), VAR, MES))
pred_art$data <- as.Date(pred_art$data)

for (nm in names(malhas)) {
  arq <- file.path(DIR, sprintf("malha_%s.rds", gsub("/", "-", nm)))
  if (file.exists(arq)) {
    if (REFAZER_COLAPSADAS && !identical(readRDS(arq)$status, "ok")) {
      file.rename(arq, sub("\\.rds$", "_colapsada.rds", arq))
      cat(sprintf("\n[%s] ajuste anterior colapsou: refazendo\n", nm))
    } else { cat(sprintf("\n[%s] ja existe, pulando\n", nm)); next }
  }
  cat(sprintf("\n===== malha %s =====\n", nm))
  mesh <- malhas[[nm]]()
  t0 <- Sys.time()
  fit <- tryCatch(bruce_fit(d, mesh, spec, window = janela,
                            theta_init = if (PARTIDA_QUENTE) novo$fit$inla$mode$theta), error = function(e) { message("ERRO: ", conditionMessage(e)); NULL })
  if (is.null(fit)) next
  minutos <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  pred <- bruce_predict(fit, grid = grid500, cap = FALSE, resample = rs)
  saveRDS(list(nome = nm, n = mesh$n, minutos = minutos, mesh = mesh,
               hyper = fit$inla$summary.hyperpar, fixed = fit$inla$summary.fixed,
               mesh_resolution = fit$mesh_resolution, status = fit$status,
               pred_mean = wrap(pred$mean)), arq)
  cat(sprintf("[%s] %d vertices, %.1f min, range %.0f m\n", nm, mesh$n, minutos,
              fit$inla$summary.hyperpar["Range for s", "mean"]))
}

# ==================================================================
# Summary
# ==================================================================
res <- lapply(list.files(DIR, "^malha_[^_]*\\.rds$", full.names = TRUE), readRDS)   # skips *_colapsada.rds
res <- res[order(vapply(res, `[[`, 0, "n"))]
fina <- unwrap(res[[length(res)]]$pred_mean)

comparar <- function(r, ref) {
  dias <- sort(unique(pred_art$data))
  dias <- dias[as.character(dias) %in% names(unwrap(r$pred_mean))]
  cmp <- vapply(dias, function(dia) {
    p <- pred_art[pred_art$data == dia, ]
    a <- terra::extract(unwrap(r$pred_mean)[[as.character(dia)]], cbind(p$x, p$y))[, 1]
    b <- terra::extract(ref[[as.character(dia)]], cbind(p$x, p$y))[, 1]
    c(cor_art = cor(a, p$pred, use = "complete"), rmse_art = sqrt(mean((a - p$pred)^2, na.rm = TRUE)),
      rmse_fina = sqrt(mean((a - b)^2, na.rm = TRUE)))
  }, numeric(3))
  rowMeans(cmp)
}
tab <- do.call(rbind, lapply(res, function(r) {
  h <- r$hyper; m <- comparar(r, fina)
  data.frame(malha = r$nome, status = r$status, vertices = r$n, minutos = round(r$minutos, 1),
             range = h["Range for s", "mean"], range_q025 = h["Range for s", "0.025quant"],
             range_q975 = h["Range for s", "0.975quant"],
             stdev_campo = h["Stdev for s", "mean"], grouprho = h["GroupRho for s", "mean"],
             aresta_tipica = r$mesh_resolution$edge, range_sobre_aresta = r$mesh_resolution$ratio,
             elev = r$fixed["elev", "mean"], cor_mapa_artigo = m["cor_art"],
             rmse_mapa_artigo = m["rmse_art"], rmse_mapa_vs_mais_fina = m["rmse_fina"])
}))
rownames(tab) <- NULL
cat("\n=== Sensibilidade a malha ===\n"); print(tab, digits = 3)
write.csv(tab, file.path(DIR, "resumo_sensibilidade_malha.csv"), row.names = FALSE)

png(file.path(DIR, "range_vs_malha.png"), width = 900, height = 600, res = 110)
par(mar = c(4.5, 4.5, 1, 1))
plot(tab$vertices, tab$range / 1000, log = "x", pch = 19, ylim = range(c(tab$range_q025, tab$range_q975)) / 1000,
     xlab = "Vertices da malha (escala log)", ylab = "Range posterior (km), media e IC 95%")
arrows(tab$vertices, tab$range_q025 / 1000, tab$vertices, tab$range_q975 / 1000, angle = 90, code = 3, length = 0.04)
text(tab$vertices, tab$range_q975 / 1000, tab$malha, pos = 3, cex = 0.8)
dev.off()
cat("Figura:", file.path(DIR, "range_vs_malha.png"), "\n")
