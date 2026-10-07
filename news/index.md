# Changelog

## bruce 0.1.0.10 (desenvolvimento)

Ferramentas para o uso no dia a dia, a partir dos testes de ponto de
orvalho e de otimização no Rio (`analises/rio/10` e `12` a `19`). Nada
muda no modelo.

- [`bruce_qc()`](https://asgouveiaa.github.io/bruce/reference/bruce_qc.md):
  controle de qualidade das estações. Compara cada estação com o fundo
  (célula mais próxima, mesmo dia), ou com a média diária das estações
  de referência quando não há fundo, e aponta as que fogem da própria
  rede (diferença média longe da mediana da rede, oscilação grande de um
  dia para o outro) ou têm poucos dias. Estações de redes de referência
  são apontadas, mas nunca retiradas. Com `drop = TRUE`, devolve os
  dados sem as estações apontadas e o relatório em `d$qc`. No Rio, teria
  pego sozinho um sensor de umidade quebrado que inflava o efeito da
  elevação em 2,5 vezes.
- [`bruce_update()`](https://asgouveiaa.github.io/bruce/reference/bruce_update.md):
  mapa de um dia novo em cerca de um minuto, com os hiperparâmetros de
  um ajuste completo anterior (aprender uma vez, mapear todo dia).
  Reajusta só uma janela curta terminando no dia, com θ fixo; campo,
  vieses, efeitos e nível do dia continuam estimados. Se o dia não tem
  dados, o mapa é uma previsão de um dia. Avisa quando o ajuste completo
  está velho (`max_age`, 45 dias por padrão): reajuste o modelo completo
  regularmente.
- Gráficos prontos:
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) de uma
  predição (um mapa por dia, média ou incerteza; `scale = "common"` para
  todos os dias na mesma escala, `"free"` para cada dia na sua própria
  escala e `"anomaly"` para cada dia menos a média da cidade naquele
  dia), [`plot()`](https://rdrr.io/r/graphics/plot.default.html) de uma
  validação cruzada (medido × previsto com intervalos e erro por estação
  ou região) e
  [`bruce_plot_series()`](https://asgouveiaa.github.io/bruce/reference/bruce_plot_series.md)
  (a cidade no mapa contra a média das estações de referência, dia a
  dia).

## bruce 0.1.0.9 (desenvolvimento)

Facilidade de uso e rigor, a partir dos testes de custo no Rio
(`analises/rio/07` a `10`):

- `bruce_cv(fit)`: a validação cruzada aceita o ajuste pronto e
  reaproveita dados, malha, especificação e janela. Cada dobra começa a
  busca do θ pelo modo do ajuste (mesmo resultado que começar do zero,
  mais rápido). Se o ajuste terminou com o GroupRho fixado após colapso,
  as dobras mantêm o mesmo valor.
- **Mudança de comportamento:**
  [`bruce_cv()`](https://asgouveiaa.github.io/bruce/reference/bruce_cv.md)
  avalia, por padrão, a leitura de uma nova estação
  (`scale = "observation"`): campo + viés da rede + efeito de uma
  estação nova + erro de medida. É a distribuição de onde vêm as
  leituras deixadas de fora, então a cobertura perto do nível nominal
  indica intervalos calibrados. A cobertura e o CRPS só do campo
  continuam no resumo (`coverage_field`, `crps_field`), e
  `scale = "field"` reproduz a avaliação do artigo. No Rio (jan/2025), a
  cobertura de 90% foi de ~0,51 no campo e ~0,85 na leitura.
- `theta_init` aceita um `bruce_fit` (em
  [`bruce_fit()`](https://asgouveiaa.github.io/bruce/reference/bruce_fit.md)
  e
  [`bruce_cv()`](https://asgouveiaa.github.io/bruce/reference/bruce_cv.md)):
  a estrutura dos dois modelos é conferida (redes com ruído estimado,
  termos, GroupRho fixo) e uma mensagem explica o que não bate. Com um
  vetor, um erro do INLA ou um tamanho diferente vira uma mensagem
  clara.
- [`bruce_check()`](https://asgouveiaa.github.io/bruce/reference/bruce_check.md)
  mostra como o ajuste foi feito: `int_strategy`, tempo, hiperparâmetros
  estimados, com partida quente (e de onde) ou fixos (e de onde), e por
  que o GroupRho foi fixado.
- Numeração de desenvolvimento (0.1.0.x) para distinguir versões
  instaladas.

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
- Opções para baratear o ajuste, todas desligadas por padrão:
  `bruce_spec(int_strategy = "eb", strategy = "gaussian")`,
  `bruce_fit(theta_init, fix_theta = TRUE)` (hiperparâmetros fixos) e,
  em
  [`bruce_cv()`](https://asgouveiaa.github.io/bruce/reference/bruce_cv.md),
  `fix_theta` e `fold_ids` (rodar só algumas dobras). O tempo do ajuste
  fica em `fit$minutes`.
- A janela do ajuste pode terminar (ou ter dias) sem nenhum dado: o
  campo espacial agora cobre todos os dias da janela (antes, um último
  dia vazio ficava de fora). Permite prever o dia seguinte.
- [`bruce_grid()`](https://asgouveiaa.github.io/bruce/reference/bruce_grid.md):
  grade de predição na resolução escolhida.
- Desvio e quantis das amostras com `matrixStats` quando instalado
  (mesmo resultado, bem mais rápido em grades finas).
- Validação:
  [`bruce_cv()`](https://asgouveiaa.github.io/bruce/reference/bruce_cv.md)
  (LSO e LRO, CRPS, cobertura, checkpoints) e
  [`bruce_simulate()`](https://asgouveiaa.github.io/bruce/reference/bruce_simulate.md)
  (cidade sintética).

Registro da validação com o Rio de Janeiro (mediana, janeiro de 2025) em
`analises/rio/README.md`.
