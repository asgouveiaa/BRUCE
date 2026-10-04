# Changelog

## bruce 0.1.0 (em desenvolvimento, repositório privado)

Primeira versão. Funções:

- Contrato de dados:
  [`bruce_data()`](https://asgouveiaa.github.io/bruce/reference/bruce_data.md)
  (observações, estações, redes de referência e enviesadas, background
  opcional) e
  [`bruce_covariates()`](https://asgouveiaa.github.io/bruce/reference/bruce_covariates.md)
  (cada covariável lida no próprio raster, na sua projeção e resolução;
  `na_value` para áreas sem dado).
- Modelo:
  [`bruce_mesh()`](https://asgouveiaa.github.io/bruce/reference/bruce_mesh.md)
  (contorno com o mesmo algoritmo do INLA),
  [`bruce_spec()`](https://asgouveiaa.github.io/bruce/reference/bruce_spec.md)
  (priors, termos, tratamento de colapso),
  [`bruce_fit()`](https://asgouveiaa.github.io/bruce/reference/bruce_fit.md)
  (detecção do colapso do campo e aviso quando a malha é grossa demais
  para o range estimado) e
  [`bruce_check()`](https://asgouveiaa.github.io/bruce/reference/bruce_check.md).
- Uso do ajuste:
  [`bruce_predict()`](https://asgouveiaa.github.io/bruce/reference/bruce_predict.md),
  [`bruce_decompose()`](https://asgouveiaa.github.io/bruce/reference/bruce_decompose.md)
  com [`plot()`](https://rdrr.io/r/graphics/plot.default.html),
  [`bruce_zonal()`](https://asgouveiaa.github.io/bruce/reference/bruce_zonal.md).
- Saídas em tabela são tibbles
  ([`bruce_cv()`](https://asgouveiaa.github.io/bruce/reference/bruce_cv.md),
  [`bruce_zonal()`](https://asgouveiaa.github.io/bruce/reference/bruce_zonal.md),
  `bruce_decompose()$variance`);
  [`as_tibble()`](https://tibble.tidyverse.org/reference/as_tibble.html)
  converte predição e decomposição em tabela longa (célula × dia) para
  `dplyr`/`ggplot2`.
- Validação:
  [`bruce_cv()`](https://asgouveiaa.github.io/bruce/reference/bruce_cv.md)
  (LSO e LRO, CRPS, cobertura, checkpoints) e
  [`bruce_simulate()`](https://asgouveiaa.github.io/bruce/reference/bruce_simulate.md)
  (cidade sintética).

Registro da validação com o Rio de Janeiro (mediana, janeiro de 2025) em
`analises/rio/README.md`.
