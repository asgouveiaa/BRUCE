# Análises

Scripts que usam o bruce com dados reais. Ficam fora do pacote (não entram no
`R CMD build`), e as saídas (`*/resultados/`) ficam fora do git.

| pasta | o que é |
|---|---|
| `rio/` | Rio de Janeiro: ajuste mensal com as redes do artigo SPDE e a validação contra o artigo |
| `berlin/` | Berlim: DWD + FAIRUrbTemp + ERA5-Land, com covariáveis globais; primeiro teste fora do Rio |

Cada script é rodado de dentro da sua pasta, por exemplo:

```bash
cd ~/bruce/analises/rio && Rscript 01_ajuste_mes.R > log_01.log 2>&1
```
