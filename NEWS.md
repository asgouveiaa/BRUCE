# bruce 0.1.0 (em desenvolvimento, repositório privado)

Primeira versão. Funções:

* Contrato de dados: `bruce_data()` (observações, estações, redes de
  referência e enviesadas, background opcional) e `bruce_covariates()`
  (cada covariável lida no próprio raster, na sua projeção e resolução;
  `na_value` para áreas sem dado).
* Modelo: `bruce_mesh()` (contorno com o mesmo algoritmo do INLA),
  `bruce_spec()` (priors, termos, tratamento de colapso), `bruce_fit()`
  (detecção do colapso do campo e aviso quando a malha é grossa demais para
  o range estimado) e `bruce_check()`.
* Uso do ajuste: `bruce_predict()`, `bruce_decompose()` com `plot()`,
  `bruce_zonal()`.
* Saídas em tabela são tibbles (`bruce_cv()`, `bruce_zonal()`,
  `bruce_decompose()$variance`); `as_tibble()` converte predição e
  decomposição em tabela longa (célula × dia) para `dplyr`/`ggplot2`.
* Opções para baratear o ajuste, todas desligadas por padrão:
  `bruce_spec(int_strategy = "eb", strategy = "gaussian")`,
  `bruce_fit(theta_init, fix_theta = TRUE)` (hiperparâmetros fixos) e, em
  `bruce_cv()`, `fix_theta` e `fold_ids` (rodar só algumas dobras). O tempo
  do ajuste fica em `fit$minutes`.
* `bruce_grid()`: grade de predição na resolução escolhida.
* Desvio e quantis das amostras com `matrixStats` quando instalado (mesmo
  resultado, bem mais rápido em grades finas).
* Validação: `bruce_cv()` (LSO e LRO, CRPS, cobertura, checkpoints) e
  `bruce_simulate()` (cidade sintética).

Registro da validação com o Rio de Janeiro (mediana, janeiro de 2025) em
`analises/rio/README.md`.
