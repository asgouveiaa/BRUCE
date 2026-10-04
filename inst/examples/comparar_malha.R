# ==================================================================
# comparar_malha.R
# Rebuilds the mesh from the bruce fit's own points with the installed bruce
# version and compares it with the article's mesh and with a mesh built by
# INLA directly from the same points. Takes seconds.
# Run from ~/scratch/spde_paper after installing bruce >= 0.1.2.
# ==================================================================

library(bruce)
library(INLA)

VAR <- "median"; MES <- "2025-01"
orig <- readRDS(sprintf("resultados_9meses_median_max/pacote_mensal_independente_%s_%s.rds", VAR, MES))
novo <- readRDS(sprintf("resultados_bruce/bruce_%s_%s.rds", VAR, MES))
d    <- novo$fit$data

cat("bruce versao:", as.character(packageVersion("bruce")),
    "| fmesher:", as.character(packageVersion("fmesher")),
    "| INLA:", as.character(packageVersion("INLA")), "\n\n")

mesh_bruce <- bruce_mesh(d, max_edge = c(5000, 40000), cutoff = 3000, offset = c(1000, 50000))

# INLA directly, as in modelo_9meses_median_max.R, from the same points
loc <- unique(rbind(unique(as.matrix(d$obs[, c("x", "y")])),
                    unique(as.matrix(d$background[, c("x", "y")]))))
bnd <- inla.nonconvex.hull(loc, convex = -0.1)
mesh_inla <- inla.mesh.2d(loc = loc, boundary = bnd, max.edge = c(5000, 40000),
                          cutoff = 3000, offset = c(1000, 50000))

cat(sprintf("\nVertices: artigo %d | INLA direto %d | bruce %d (antes do ajuste: %d)\n",
            orig$mesh$n, mesh_inla$n, mesh_bruce$n, novo$fit$mesh$n))
cat("Malha do bruce identica a do artigo:",
    isTRUE(all.equal(unname(mesh_bruce$loc[, 1:2]), unname(orig$mesh$loc[, 1:2]))), "\n")
