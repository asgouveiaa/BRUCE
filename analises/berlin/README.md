# Berlim

Rede de referência: DWD (baixada automaticamente dos dados abertos).
Rede enviesada: FAIRUrbTemp Berlim (TU Berlin, 11 estações, jan/2020–jan/2024;
baixar de https://doi.org/10.48620/93247 e apontar `DIR_FAIR`).
Background: ERA5-Land horário (precisa da chave do CDS no `ecmwfr`).
Covariáveis: elevação (elevatr), fração construída (ESA WorldCover), SVF (LCZ4r).

```bash
cd ~/bruce/analises/berlin && Rscript 01_ajuste_mes.R > log_berlim.log 2>&1
```
