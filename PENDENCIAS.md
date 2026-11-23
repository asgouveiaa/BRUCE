# Pendências antes de tornar o repositório público

## Validação
- [ ] Teste de sensibilidade à malha no Rio (`analises/rio/06_sensibilidade_malha.R`):
      o range converge com malhas mais finas?
- [ ] Rever a heurística padrão do `bruce_mesh()` à luz desse teste (hoje gera
      malhas grossas para redes com o espaçamento da do Rio).
- [ ] Conferir em ajuste real o aviso de resolução da malha (`fit$mesh_resolution`);
      só foi testado com uma malha de exemplo.
- [ ] LSO e LRO na mínima do Rio.
- [ ] Berlim (`analises/berlin/01_ajuste_mes.R`): primeiro teste fora do Rio.
- [ ] Um segundo caso fora do Rio (Zurique ou Birmingham).

## Custo
- [x] Rodar `analises/rio/07_custo.R` (depois do teste da malha) e decidir:
      solver/threads padrão, CV com partida quente ou θ fixo, mapa diário
      operacional com θ fixo.
- [ ] Resolução da grade (500 m a 5 km): predizer direto na grade grossa ou
      predizer a 500 m e agregar? Decidir com a parte C do 07.
- [ ] Medir o pico de RAM (`/usr/bin/time -v`) por jeito de ajustar (ccd, eb,
      θ fixo) e por malha; tamanho do objeto do ajuste com `keep_config`.
      Documentar os requisitos de máquina.
- [x] Testar θ emprestado entre meses (`08`, `10`): nov/2024, jan/2024 e jun/2025 →
      jan/2025; mapa, erro e cobertura contra o ajuste completo.
- [ ] Repetir a parte B do 07 em outro mês (inverno; máxima e mínima).
- [x] `bruce_cv(fit)` com partida quente; avaliação na escala da leitura
      (padrão) e do campo; trava e procedência do θ emprestado; `bruce_check`
      transparente (0.1.0.9).
- [ ] `eb` como padrão do `bruce_fit` (antes: 08 parte 1 com a máxima).
- [x] Função do mapa diário com θ emprestado: `bruce_update()` (0.1.0.10).
- [ ] `eb` com partida quente como modo rápido do `bruce_cv` (testado só em
      jan/2025, 6 estações).
- [ ] Mudanças no modelo para baratear (depois do 07): malha não uniforme,
      fixar sd de redes estáveis, janelas mais curtas.

## Testes automáticos
- [ ] Teste de `bruce_cv(scheme = "lro")` (hoje só o LSO é testado).
- [ ] Teste de recuperação dos parâmetros com mais dias e redes aglomeradas
      (`bruce_simulate(clustered = TRUE)`), incluindo o caminho de colapso.
- [ ] Teste de `bruce_predict()` com rasters em projeções diferentes e `resample = "near"`.
- [ ] Teste de `bruce_zonal()` com pesos.
- [ ] `R CMD check` limpo e um fluxo de CI que rode os testes sem INLA.

## Documentação
- [ ] Gráficos prontos no tutorial: `plot()` da predição e da validação,
      `bruce_plot_series()` (0.1.0.10).
- [ ] Figuras de um ajuste real no tutorial (Rio ou Berlim).
- [ ] Função para CPO/PIT (usada no resumo do congresso).
- [ ] Função para rodar várias janelas com checkpoint (os 9 meses do artigo).
- [ ] CITATION e DOI (Zenodo) quando abrir.

## Mudanças no pacote que vieram dos testes de ponto de orvalho (06/10/2026)
Nada disso foi implementado ainda; os testes estão em `analises/rio/12` a `16`.
- [x] Controle de qualidade de estações contra o fundo: `bruce_qc()`
      (0.1.0.10; antes `qc=3` no script 12): tirar estações de redes cidadãs cuja diferença ao ERA5 fica
      longe da diferença típica da própria rede, que oscilam demais, ou com
      poucos dias. Nunca tirar estações de referência. Um sensor quebrado
      (wunder:IRIODE84) inflava sozinho o efeito da elevação em 2,5 vezes.
- [ ] Buffer do fundo: avisar quando as células do fundo superam muito as
      observações de estação ou vão longe da área das estações (com 40 km, o
      ERA5 dominava o ponto de orvalho e criava um range de ~400 km; com 10 km,
      tudo melhorou). Pensar num padrão ligado ao tamanho da área das estações.
- [ ] Viés próprio para o fundo (opcional, desligado por padrão): o ERA5-Land
      d2m é 0,2–0,9 °C mais seco que INMET e aeroportos em 2023 e 2026.
- [ ] sd do efeito de estação por tipo de rede (referência x cidadã): com um
      sd só, a cobertura na escala da leitura fica em ~1,00.
- [ ] Valor de correlação entre dias para o colapso: emprestar de um mês com
      muitos dados, em vez do 0,55 genérico (no ponto de orvalho, 0,9 veio
      de jan/2026).
- [ ] Função de índice de calor a partir de dois ajustes (T e Td), com
      incerteza por sorteios (protótipo: `analises/rio/16_indice_calor.R`).
- [ ] Documentar: ponto de orvalho é experimental. Com poucas estações (2023:
      19), o campo só capta escala regional (range ~160 km, não identificado);
      o detalhe local vem das covariáveis.

## Ideias futuras (depois da primeira versão pública)
- Umidade: modelar o ponto de orvalho (Td, em °C), não a umidade relativa
  (limitada a 0–100% e dependente de T). ERA5-Land tem `d2m` para o fundo.
  Covariáveis candidatas: elevação, distância ao mar/água, vegetação.
  Primeiro levantar quais redes medem umidade e a qualidade (sensores baratos
  de UR têm viés e saturam em 100%).
- Índice de calor a partir dos mapas de T e Td (UR por Magnus; índice do NWS
  ou Humidex). Cuidados: com dados diários, usar uma aproximação documentada
  (ex.: índice(Tmax, Td média)); incerteza por amostras de T e Td, que, se
  tratadas como independentes, subestima a do índice. Solução completa:
  modelo conjunto T–Td.
- Previsão: (1) persistência pelo AR1 do próprio modelo (teste: `09`);
  (2) downscaling de uma previsão do tempo em grade usada como fundo, com
  correção de viés; meio-termo: nível da cidade da previsão oficial +
  padrão espacial do bruce.
