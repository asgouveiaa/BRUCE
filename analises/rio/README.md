# Rio de Janeiro

Antes de rodar: ajuste os caminhos em `config.R` e ponha a chave do
Wunderground em `~/.Renviron` (`WU_API_KEY=...`). Rode tudo de dentro desta pasta.

| script | o que faz | tempo |
|---|---|---|
| `query_meteorj_periodo.R` | consulta às redes (INMET, tb_clima_ar, REDEMET, Aurassure, Wunderground); `estatisticas = c("median", "max", "min")` | — |
| `01_ajuste_mes.R` | ajusta um mês (`VAR` = median, max ou min), prediz a 500 m, decompõe; LRO opcional | ~25 min |
| `02_comparar_com_artigo.R` | compara o ajuste do bruce com o salvo do artigo (efeitos, hiperparâmetros, mapas) | segundos |
| `03_comparar_entradas.R` | compara as entradas (estações, valores, covariáveis, ERA5, malha) sem reajustar | segundos |
| `04_comparar_entradas_identicas.R` | ajusta o bruce com as entradas exatas do `.rds` do artigo | ~20 min |
| `05_comparar_malha.R` | refaz a malha e compara com a do artigo | segundos |
| `06_sensibilidade_malha.R` | mesmo modelo com 4–5 malhas, da mais grossa à mais fina | horas |

## Registro da validação (mediana, janeiro de 2025)

- Com as entradas e a malha do artigo (`04`), o bruce reproduz o ajuste:
  efeitos fixos a menos de 0,03 desvio posterior; range 7.682 m contra 7.669 m.
- Com as próprias entradas (`01`) e a mesma malha do artigo: range 7.521 m.
- Com a malha refeita pelo fmesher atual (1.022 vértices contra 1.013 do artigo,
  mesmos parâmetros): range 4.146 m, mas mapas com correlação 0,99 e RMSE 0,16 °C
  em relação ao artigo. O range é sensível à malha; ver `06`.
- Mínima (janeiro de 2025): viés do Wunderground ~0, SVF com efeito
  (−0,41; IC −0,75 a −0,07), range 2,3 km (abaixo da aresta da malha).
