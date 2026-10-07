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
| `06b_resumo_malha.R` | tabela e figura do 06 com as malhas já terminadas (não reajusta) | segundos |
| `07_custo.R` | custo: formas de rodar o INLA (threads, PARDISO, `eb`, θ fixo), CV barata e resolução da grade (1000 a 100 m) | ~8 h |
| `08_theta_emprestado.R` | as janelas do artigo com θ próprio (`eb`) e com θ emprestado de cada outra janela: os achados se mantêm? | horas |
| `09_previsao_dia_seguinte.R` | previsão do dia seguinte a cada janela do artigo (persistência do modelo) contra "amanhã = hoje" e a média do mês | ~15–30 min |
| `09b_mapas_previsao.R` | mapas do último dia, da previsão e da incerteza para cada janela do 09 | ~15 min |
| `10_cobertura_theta_emprestado.R` | LSO com θ próprio (reaprendido), fixo e emprestado de outras janelas: erro e cobertura dos intervalos | ~2,5 h |
| `11_levantamento_umidade.R` | quais redes medem umidade/ponto de orvalho e a qualidade dos dados | segundos |
| `query_orvalho_periodo.R` | ponto de orvalho diário por estação (INMET, aeroportos, Alerta Rio, Aurassure, Wunderground) | — |
| `12_orvalho_ajuste.R` | primeiro ajuste do ponto de orvalho com o bruce (um mês) + LSO | ~30–60 min |

## Registro da validação (mediana, janeiro de 2025)

- Com as entradas e a malha do artigo (`04`), o bruce reproduz o ajuste:
  efeitos fixos a menos de 0,03 desvio posterior; range 7.682 m contra 7.669 m.
- Com as próprias entradas (`01`) e a mesma malha do artigo: range 7.521 m.
- Com a malha refeita pelo fmesher atual (1.022 vértices contra 1.013 do artigo,
  mesmos parâmetros): range 4.146 m, mas mapas com correlação 0,99 e RMSE 0,16 °C
  em relação ao artigo. O range é sensível à malha; ver `06`.
- Mínima (janeiro de 2025): viés do Wunderground ~0, SVF com efeito
  (−0,41; IC −0,75 a −0,07), range 2,3 km (abaixo da aresta da malha).
