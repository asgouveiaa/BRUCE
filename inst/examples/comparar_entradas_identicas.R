# ==================================================================
# comparar_entradas_identicas.R
# Fits bruce on EXACTLY the inputs saved in the article's package for one
# month (stations, ERA5 cells, standardised covariates and mesh), so the
# only thing that differs from the article is the modelling code.
#
#   Range close to the article's  -> bruce reproduces the model; earlier
#                                    differences come from the inputs
#                                    (covariate rebuilding in rio_adapter.R).
#   Range still different         -> a real difference in the model code.
#
# Run from ~/scratch/spde_paper (about the same time as one month).
# ==================================================================

library(bruce)
library(sf)

VAR <- "median"; MES <- "2025-01"
orig  <- readRDS(sprintf("resultados_9meses_median_max/pacote_mensal_independente_%s_%s.rds", VAR, MES))
fit_o <- orig$fit

# ---- 0. Two cheap checks first ---------------------------------------------
cat("Range posterior (artigo):\n")
print(fit_o$summary.hyperpar["Range for s", c("mean", "sd", "0.025quant", "0.975quant")])
novo <- readRDS(sprintf("resultados_bruce/bruce_%s_%s.rds", VAR, MES))
cat("Range posterior (bruce via rio_adapter):\n")
print(novo$fit$inla$summary.hyperpar["Range for s", c("mean", "sd", "0.025quant", "0.975quant")])
cat(sprintf("Vertices da malha: artigo %d | bruce %d\n\n", orig$mesh$n, novo$fit$mesh$n))

# ---- 1. Article's data, mapped onto the bruce contract -----------------------
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
                      bias   = c(FALSE, TRUE, TRUE, TRUE),
                      obs_sd = c(0.3, 0.3, NA, NA))

d <- bruce_data(obs, stations, sources, crs = 31983, background = background, background_sd = 0.7)

# ---- 2. The article's standardised covariates, injected as they are ----------
st_cov <- data.frame(station_id = est$estacao_id, elev = est$elev_std, urb = est$urb_std,
                     svf = est$svf_std)
st_cov <- st_cov[!duplicated(st_cov$station_id), ]
st_cov$svf[is.na(st_cov$svf)] <- 0            # INLA treats a missing covariate as 0
bg_cov <- data.frame(x = bg$x, y = bg$y, elev = bg$elev_std, urb = bg$urb_std, svf = 0)
bg_cov <- bg_cov[!duplicated(bg_cov[, c("x", "y")]), ]
d$covariates <- list(names = c("elev", "urb", "svf"), rasters = NULL, scaling = NULL,
                     station = st_cov, background = bg_cov)

# ---- 3. Same mesh, same priors -----------------------------------------------
spec <- bruce_spec(range_prior = 15000, intercept_prior = c(0, 0.001), on_collapse = "warn")
t0 <- Sys.time()
fit_x <- bruce_fit(d, orig$mesh, spec)
message(sprintf("Ajuste: %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# ---- 4. Compare ---------------------------------------------------------------
hp <- function(f, p) { i <- grep(p, rownames(f$summary.hyperpar)); if (length(i)) f$summary.hyperpar[i[1], "mean"] else NA }
cmp <- data.frame(
  parametro = c("Range (m)", "Stdev campo", "GroupRho", "Rho diario", "Prec estacao"),
  artigo = c(hp(fit_o, "^Range for s"), hp(fit_o, "^Stdev for s"), hp(fit_o, "^GroupRho for s"),
             hp(fit_o, "^Rho for dia_idx"), hp(fit_o, "^Precision for id_estacao_idx")),
  bruce_mesmas_entradas = c(hp(fit_x$inla, "^Range for s"), hp(fit_x$inla, "^Stdev for s"),
                            hp(fit_x$inla, "^GroupRho for s"), hp(fit_x$inla, "^Rho for day_idx"),
                            hp(fit_x$inla, "^Precision for station_idx")),
  bruce_adapter = c(hp(novo$fit$inla, "^Range for s"), hp(novo$fit$inla, "^Stdev for s"),
                    hp(novo$fit$inla, "^GroupRho for s"), hp(novo$fit$inla, "^Rho for day_idx"),
                    hp(novo$fit$inla, "^Precision for station_idx")))
cat("\n=== Hiperparametros ===\n"); print(cmp, digits = 4, row.names = FALSE)

nm <- c(b0 = "intercept", elev_std = "elev", urb_std = "urb", svf_std = "svf",
        aura = "bias_aurassure", wunder = "bias_wunder", aeroporto = "bias_airport")
nm <- nm[names(nm) %in% rownames(fit_o$summary.fixed) & nm %in% rownames(fit_x$inla$summary.fixed)]
fx <- data.frame(termo = names(nm),
                 artigo = fit_o$summary.fixed[names(nm), "mean"],
                 bruce_mesmas_entradas = fit_x$inla$summary.fixed[nm, "mean"],
                 diff_em_sd = (fit_x$inla$summary.fixed[nm, "mean"] - fit_o$summary.fixed[names(nm), "mean"]) /
                   fit_o$summary.fixed[names(nm), "sd"])
cat("\n=== Efeitos fixos ===\n"); print(fx, digits = 3, row.names = FALSE)
cat("\nStatus:", fit_x$status, "\n")
