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
* Validação: `bruce_cv()` (LSO e LRO, CRPS, cobertura, checkpoints) e
  `bruce_simulate()` (cidade sintética).

Registro da validação com o Rio de Janeiro (mediana, janeiro de 2025) em
`analises/rio/README.md`.
