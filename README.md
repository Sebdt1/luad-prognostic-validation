# luad-prognostic-validation

¿Un score transcriptómico de proliferación aporta información pronóstica más allá
del estadio, la edad y el sexo en adenocarcinoma de pulmón? ¿Y replica en una
cohorte independiente?

Aquí está el estudio entero: pipeline reproducible en R/Bioconductor sobre
TCGA-LUAD y GSE68465, con validación interna y externa. El reporte no infla los
resultados, que es buena parte de la gracia.

## Qué salió

El score sí se asocia de forma independiente a la supervivencia tras ajustar por
estadio, edad y sexo (HR 1.35, IC95% 1.11–1.63, p=0.003), y la asociación replica
en GEO (HR 1.24, IC95% 1.08–1.42, p=0.002). Hasta ahí, bien.

Lo que no se sostiene es la mejora de discriminación. El C-index corregido por
optimismo sube de 0.66 a 0.69, pero los contrastes formales de ΔAUC y ΔBrier en el
test retenido no son significativos: los intervalos cruzan cero. Un elastic-net de
26 genes, que en training parecía muy superior (C=0.77), se hunde a 0.56 en el test
retenido, con un intervalo que incluye 0.5. Sobreajuste de manual. Y al transportar
el modelo congelado a GEO el score empeora la calibración de forma significativa,
aunque re-estandarizarlo en la cohorte destino lo corrige del todo.

Señal real y replicable, aporte pequeño, y la complejidad no ayuda. El valor del
trabajo está en la validación, no en la firma.

Los resultados completos, con tablas y 9 figuras, están en
[`results/review_packet.md`](results/review_packet.md), que GitHub renderiza
directamente.

## Estado

Corrida completa ejecutada el 2026-07-18 sobre TCGA-LUAD (n=504, 182 muertes) y
GSE68465 (n=442, 236 muertes). El pipeline corre sin error de principio a fin, los
22 tests de `tests/test_utils.R` pasan y `renv.lock` está sincronizado. Entorno
verificado: R 4.6.1, Bioconductor 3.23, Quarto 1.9.38.

### Cosas que costaron y conviene saber

La descarga de RNA-seq de TCGAbiolinks 2.40.0 está rota por sus dos vías. Con
`method="api"`, la función interna `GDCdownload.aux()` termina en `message()`, que
devuelve `NULL`, y el llamador hace `if (ret == 1)`: error de argumento de longitud
cero. Con `method="client"` descarga el zip de gdc-client pero no lo extrae ni lo
ejecuta. La solución que quedó en `R/02_acquire.R` usa `GDCquery` solo para el
catálogo, descarga con gdc-client vía manifiesto y construye el
`SummarizedExperiment` leyendo los `.tsv` STAR. Los datos son los mismos del GDC.

Tres detalles más del entorno, por si te ahorran una tarde:

- El timeout de descarga por defecto de R son 60 segundos y trunca los archivos
  grandes en torno a los 250 MB. Aquí se sube a 3600.
- La descarga va a un directorio local (`%LOCALAPPDATA%/luad_gdc`) y no al
  proyecto, porque OneDrive interfiere al sincronizar miles de archivos.
- Si `quarto` no está en el PATH, `tar_quarto` falla al cargar el pipeline; basta
  con exportar `QUARTO_PATH`. Y si la librería de R por defecto no es escribible,
  crea la personal antes de correr `bootstrap.R`.

---

## Cómo reproducir

```bash
# 1) Entorno reproducible (instala stack CRAN + Bioconductor y crea renv.lock)
Rscript setup/bootstrap.R

# 2) (En corridas posteriores, en otra máquina) restaurar el entorno exacto
Rscript -e 'renv::restore()'

# 3) Correr el pipeline (respeta run.smoke_test de config/params.yml)
Rscript -e 'targets::tar_make()'

# 4) Ver el DAG / estado
Rscript -e 'targets::tar_visnetwork()'
Rscript -e 'targets::tar_progress()'

# 5) Tests unitarios de las utilidades críticas
Rscript -e 'testthat::test_file("tests/test_utils.R")'
```

El **reporte científico** (`reports/report.qmd`) se renderiza como parte del
pipeline (target `report`, vía `tarchetypes::tar_quarto`). También manual:

```bash
quarto render reports/report.qmd
```

### Modo smoke vs. completo

Todo se controla en `config/params.yml`:

- `run.smoke_test: true` → submuestrea la cohorte (estratificado y determinista)
  para iterar rápido. La **primera descarga** de TCGA/GEO es inevitable pero se
  **cachea** en `data/raw/` (idempotente); el smoke acelera el cómputo posterior.
- `run.smoke_test: false` → corrida completa.

---

## Estructura

```
luad-prognostic-validation/
├── README.md                 # este archivo
├── .Rprofile                 # activa renv si ya está inicializado
├── .gitignore                # data/raw, renv/library, artefactos, etc.
├── config/
│   └── params.yml            # ÚNICA fuente de verdad de parámetros + semilla
├── setup/
│   ├── dependencies.R        # manifiesto de paquetes (CRAN + Bioconductor)
│   └── bootstrap.R           # instala el stack y GENERA renv.lock (snapshot)
├── _targets.R                # pipeline reproducible (DAG)
├── R/
│   ├── 00_utils.R            # utilidades puras (semilla, escalas, estadio, IO)
│   ├── 01_split.R            # partición train/test ANTES de todo análisis
│   ├── 02_acquire.R          # TCGA RNA-seq + TCGA-CDR + cohorte GEO (cacheado)
│   ├── 03_qc.R               # filtro de genes, VST/logCPM, QC, PCA de lote
│   ├── 04_eda.R              # cohorte, faltantes, KM global/estratificado
│   ├── 05_score.R            # singscore (G2M/E2F) + Cox elastic-net (train)
│   ├── 06_model.R            # Cox incremental + supuestos + validación interna
│   └── 07_external_valid.R   # transporte y replicación en GEO
├── data/
│   ├── raw/                  # descargas (en .gitignore)
│   └── processed/            # intermedios (en .gitignore)
├── results/
│   ├── figures/              # figuras autogeneradas
│   └── tables/               # tablas autogeneradas (CSV)
├── reports/
│   └── report.qmd            # reporte científico (Quarto)
└── tests/
    └── test_utils.R          # tests testthat de la lógica crítica
```

---

## Decisiones metodológicas clave

- **Partición antes de analizar y selección solo en training.** Evita fuga de
  información: el filtro de genes, la escala del score y el elastic-net se
  aprenden solo en *training*; el *test* queda intacto para validación interna.
- **Endpoint curado (TCGA-CDR), no `days_to_death`.** OS/PFI del recurso curado
  de Liu et al. (Cell 2018) son más fiables que los campos clínicos crudos.
- **Ajuste obligatorio por estadio.** El estadio es el confusor pronóstico
  dominante en LUAD; sin ajustarlo, una firma proliferativa suele ser su proxy.
- **Score rank-based (singscore) para transferir entre plataformas.** La
  discriminación es invariante a re-escalados monótonos, así que el score se
  transporta a GEO sin re-entrenar; la escala z solo hace interpretable el HR.
- **VST con dispersión ajustada en training.** La normalización no usa el
  desenlace (no es fuga selectiva), pero ajustar en training aporta rigor extra.
- **Lote se documenta, no se sobre-corrige.** Una corrección agresiva puede
  borrar señal biológica confundida con el lote.
- **Validación honesta.** C-index corregido por optimismo (bootstrap), test
  interno retenido, AUC/t, calibración y Brier/IPA con IC; y validación externa
  que separa *discriminación* (suele conservarse) de *calibración* (suele
  degradarse → recalibrar el hazard basal, no descartar el marcador).
- **`renv.lock` generado, no escrito a mano.** Un lockfile válido requiere
  versiones y hashes reales; se captura con `renv::snapshot()` en `bootstrap.R`.
- **RSF disponible pero no ejecutado.** `R/06_model.R::fit_rsf()` implementa el
  random survival forest como sensibilidad, pero **no se incluye como target**:
  el protocolo lo contempla solo *si falla la proporcionalidad de riesgos*, y en
  la cohorte completa el test de Schoenfeld **no la rechaza** (estadio p=0.10,
  global p=0.052). Las sensibilidades que sí corren son Cox estratificado, RMST
  e imputación múltiple (`mice`).

---

## Amenazas a la validez (resumen)

Sesgo de espectro (cohorte quirúrgica), confusión residual (tabaquismo/
tratamiento mal capturados), ancestría (TCGA sesgada a ascendencia europea),
efecto de lote/plataforma, multiplicidad y sobreajuste, faltantes (caso completo
vs. `mice`), y fiabilidad del endpoint. Se discuten en `reports/report.qmd`.

---

## Datos y licencias

- **TCGA-LUAD** (RNA-seq STAR counts) y **TCGA-CDR**: open-access del GDC (sin
  token). Sujetos a los términos de uso de datos de TCGA/GDC.
- **GEO** (`GSE68465` por defecto; alternativa `GSE31210`): datos públicos.
- Los datos crudos **no se versionan** (`data/raw/` en `.gitignore`); el pipeline
  los descarga de forma idempotente.
