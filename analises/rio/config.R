# ==================================================================
# analises/rio/config.R
# Caminhos usados por todos os scripts do Rio. Ajuste aqui, uma vez.
# Os scripts devem ser rodados de dentro de analises/rio:
#   cd ~/bruce/analises/rio && Rscript 01_ajuste_mes.R
# ==================================================================

DIR_ARTIGO   <- path.expand("~/scratch/spde_paper")                    # pasta do artigo SPDE
ARTIGO       <- file.path(DIR_ARTIGO, "resultados_9meses_median_max")  # ajustes salvos do artigo
MALHA_AP     <- file.path(DIR_ARTIGO, "malhas", "ap_malha.gpkg")       # regioes para o LRO
WU_CACHE     <- file.path(DIR_ARTIGO, "wu_cache_v2.duckdb")            # cache do Wunderground ja existente

ERA5_FUNCOES <- path.expand("~/MEGA/cisc/era5_interp/com_aura_novo/funcoes downscaling ERA5 MRJ.R")
DIR_ERA5     <- path.expand("~/MEGA/cisc/era5_interp/independente")
USO_SOLO     <- path.expand("~/MEGA/cisc/era5_interp/Uso_do_Solo_2019.shp")

RES <- "resultados"   # saidas do bruce (relativo a analises/rio; fora do git)
dir.create(RES, showWarnings = FALSE, recursive = TRUE)

# A chave da API do Wunderground fica fora do codigo, em ~/.Renviron:
#   WU_API_KEY=sua_chave
if (!nzchar(Sys.getenv("WU_API_KEY"))) {
  warning("WU_API_KEY nao definida em ~/.Renviron: dias do Wunderground fora do cache nao serao baixados")
}
