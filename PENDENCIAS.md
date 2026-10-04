# Pendências antes de tornar o repositório público

## Validação

Teste de sensibilidade à malha no Rio
(`analises/rio/06_sensibilidade_malha.R`): o range converge com malhas
mais finas?

Rever a heurística padrão do
[`bruce_mesh()`](https://asgouveiaa.github.io/bruce/reference/bruce_mesh.md)
à luz desse teste (hoje gera malhas grossas para redes com o espaçamento
da do Rio).

Conferir em ajuste real o aviso de resolução da malha
(`fit$mesh_resolution`); só foi testado com uma malha de exemplo.

LSO e LRO na mínima do Rio.

Berlim (`analises/berlin/01_ajuste_mes.R`): primeiro teste fora do Rio.

Um segundo caso fora do Rio (Zurique ou Birmingham).

## Testes automáticos

Teste de `bruce_cv(scheme = "lro")` (hoje só o LSO é testado).

Teste de recuperação dos parâmetros com mais dias e redes aglomeradas
(`bruce_simulate(clustered = TRUE)`), incluindo o caminho de colapso.

Teste de
[`bruce_predict()`](https://asgouveiaa.github.io/bruce/reference/bruce_predict.md)
com rasters em projeções diferentes e `resample = "near"`.

Teste de
[`bruce_zonal()`](https://asgouveiaa.github.io/bruce/reference/bruce_zonal.md)
com pesos.

`R CMD check` limpo e um fluxo de CI que rode os testes sem INLA.

## Documentação

Figuras de um ajuste real no tutorial (Rio ou Berlim).

Função para CPO/PIT (usada no resumo do congresso).

Função para rodar várias janelas com checkpoint (os 9 meses do artigo).

CITATION e DOI (Zenodo) quando abrir.
