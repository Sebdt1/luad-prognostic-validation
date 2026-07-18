# Review packet — validación pronóstica TCGA-LUAD + GEO

Generado: 2026-07-18 17:37 · Documento para revisión externa.
Todas las cifras provienen del store de `targets` de la corrida registrada abajo.

## 1. Árbol del proyecto (2 niveles)

```
.
├── _targets.R
├── config/
│   ├── params.yml
├── data/
│   ├── processed
│   ├── raw
├── R/
│   ├── 00_utils.R
│   ├── 01_split.R
│   ├── 02_acquire.R
│   ├── 03_qc.R
│   ├── 04_eda.R
│   ├── 05_score.R
│   ├── 06_model.R
│   ├── 07_external_valid.R
│   ├── 08_extra_analyses.R
│   ├── 09_figures.R
├── README.md
├── renv.lock
├── reports/
│   ├── report.html
│   ├── report.qmd
├── results/
│   ├── figures
│   ├── review_packet.md
│   ├── tables
├── setup/
│   ├── bootstrap.R
│   ├── dependencies.R
├── tests/
│   ├── test_utils.R
├── renv/            (librería del proyecto; no versionada)
└── _targets/        (store del pipeline; no versionado)
```

## 2. config/params.yml (completo)

```yaml
# =============================================================================
# params.yml — Configuración central y única fuente de verdad del pipeline
# -----------------------------------------------------------------------------
# Todo parámetro que afecte un resultado vive aquí (no incrustado en el código),
# para que el pipeline sea auditable y reproducible. La semilla global se fija
# desde aquí y se reporta en sessionInfo(). Cambiar 'run.smoke_test' a false
# ejecuta la corrida completa; en true se submuestrea para iterar rápido.
# =============================================================================

project:
  name: luad-prognostic-validation
  # Semilla global. Se aplica en R/00_utils.R::set_global_seed() y se reporta.
  seed: 1234

run:
  # smoke_test = true  -> submuestra las cohortes para validar el DAG en minutos.
  # smoke_test = false -> corrida completa (descargas grandes, horas la 1a vez).
  smoke_test: false
  # Nº de muestras por brazo (train/test) que se conservan en modo smoke.
  smoke_n: 80
  cache_dir: data/raw
  processed_dir: data/processed

split:
  # La partición se hace ANTES de cualquier análisis. Toda selección de
  # features/score se realiza SOLO en training para evitar fuga de información.
  p_train: 0.70
  # Estratificamos por estadio agrupado para equilibrar el pronóstico basal
  # entre brazos (el estadio es el confusor pronóstico dominante en LUAD).
  stratify_by: stage_group

acquisition:
  tcga_project: TCGA-LUAD
  gdc_workflow: "STAR - Counts"
  # TP = Primary Tumor. Excluimos normales y metástasis para una cohorte
  # pronóstica homogénea (documentado como sesgo de espectro: cohorte quirúrgica).
  sample_type_codes: ["01"]
  cdr:
    # Recurso clínico curado TCGA-CDR (Liu et al., Cell 2018). Usamos OS/PFI
    # curados, NO days_to_death crudo (menos fiable). El UUID es de la página
    # PanCanAtlas del GDC; VERIFICAR en:
    #   https://gdc.cancer.gov/about-data/publications/pancanatlas
    # Si la descarga automática falla, descargar manualmente el archivo
    # 'TCGA-CDR-SupplementalTableS1.xlsx' a data/raw/ y el pipeline lo usará.
    url: "https://api.gdc.cancer.gov/data/1b5f413e-a8d1-4d10-92eb-7c4ae739ed81"
    file: TCGA-CDR-SupplementalTableS1.xlsx
    sheet: "TCGA-CDR"
  # Endpoint pronóstico primario. OS = supervivencia global; PFI = intervalo
  # libre de progresión (recomendado por TCGA-CDR para LUAD junto con OS).
  endpoint: OS
  geo:
    # Cohorte de validación externa independiente (adyuvante/quirúrgica).
    # GSE68465 (Director's Challenge, Affymetrix U133A) trae supervivencia.
    # Alternativa: GSE31210 (estadio I-II, japonesa).
    gse_id: GSE68465
    platform: GPL96

qc:
  # Filtro de genes de baja expresión (estilo edgeR::filterByExpr, criterio
  # explícito y auditable): conservar genes con >= min_count en >= min_prop
  # de las muestras del TRAINING (el umbral se aprende solo en training).
  min_count: 10
  min_prop_samples: 0.20
  # Transformación de la varianza para análisis downstream y visualización.
  # vst (DESeq2) recomendado; logcpm (edgeR) como alternativa.
  transform: vst
  # Detección de muestras atípicas por tamaño de librería (nº de MADs).
  library_size_mad: 4
  # Variables candidatas a efecto de lote a inspeccionar por PCA (documentar,
  # NO sobre-corregir): centro de tejido (TSS) y placa.
  batch_vars: ["tss", "plate"]

score:
  # (a) Score primario de proliferación con singscore (rank-based, transferible
  #     entre plataformas sin re-entrenamiento) sobre firmas Hallmark.
  method: singscore
  gene_sets:
    - HALLMARK_G2M_CHECKPOINT
    - HALLMARK_E2F_TARGETS
  msigdb_species: "Homo sapiens"
  msigdb_collection: H
  # Escalado para transferencia entre plataformas: z-score respecto a la
  # distribución del TRAINING (los parámetros de escala se guardan y se aplican
  # tal cual a test y a GEO -> sin fuga).
  scaling: zscore
  # (b) Score secundario: Cox elastic-net (glmnet), lambda por CV.
  elasticnet:
    alpha: 0.5
    nfolds: 10
    # nº mínimo de eventos para intentar el elastic-net (si no, se omite).
    min_events: 40

model:
  # Modelo clínico estándar. SIEMPRE se ajusta por estadio: una "firma" sin
  # ajuste por estadio suele ser un mero proxy del estadio.
  clinical_vars: ["age_std", "gender", "stage_group"]
  time_var: time
  event_var: status

validation:
  # C-index corregido por optimismo vía bootstrap.
  bootstrap_B: 200
  # Horizontes (días) para AUC tiempo-dependiente (timeROC) y Brier/IPA: 1,3,5 años.
  horizons_days: [365, 1095, 1825]
  # Intervalo de confianza para todas las métricas reportadas.
  ci_level: 0.95

report:
  # Umbral de p para señalar (no para dicotomizar decisiones): solo descriptivo.
  alpha: 0.05
```

## 3. Reproducibilidad

- **Modo de corrida**: `run.smoke_test = false` → **COHORTE COMPLETA**
- **targets**: 69 targets registrados; `tar_outdated()` = vacío (pipeline al día)
- **renv**: consistent state (lockfile sincronizado)
- **Semilla global**: 1234

### sessionInfo()

```
R version 4.6.1 (2026-06-24 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 11 x64 (build 26200)

Matrix products: default
  LAPACK version 3.12.1

locale:
[1] LC_COLLATE=Spanish_Mexico.utf8  LC_CTYPE=Spanish_Mexico.utf8   
[3] LC_MONETARY=Spanish_Mexico.utf8 LC_NUMERIC=C                   
[5] LC_TIME=Spanish_Mexico.utf8    

time zone: America/Bogota
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices datasets  utils     methods   base     

other attached packages:
[1] timeROC_0.4.1             riskRegression_2026.03.11
[3] survival_3.8-9            targets_1.12.0           

loaded via a namespace (and not attached):
 [1] tidyselect_1.2.1       dplyr_1.2.1            farver_2.1.2          
 [4] S7_0.2.2               fastmap_1.2.0          TH.data_1.1-5         
 [7] digest_0.6.39          rpart_4.1.27           base64url_1.4         
[10] lifecycle_1.0.5        secretbase_1.3.0       cluster_2.1.8.2       
[13] processx_3.9.0         magrittr_2.0.5         compiler_4.6.1        
[16] rlang_1.3.0            Hmisc_5.2-6            tools_4.6.1           
[19] igraph_2.3.3           yaml_2.3.12            data.table_1.18.4     
[22] knitr_1.51             prettyunits_1.2.0      timereg_2.0.7         
[25] htmlwidgets_1.6.4      RColorBrewer_1.1-3     multcomp_1.4-31       
[28] polspline_1.1.25       withr_3.0.3            foreign_0.8-91        
[31] numDeriv_2016.8-1.1    pec_2025.06.24         nnet_7.3-20           
[34] grid_4.6.1             mets_1.3.11            colorspace_2.1-3      
[37] future_1.70.0          ggplot2_4.0.3          globals_0.19.1        
[40] scales_1.4.0           iterators_1.0.14       MASS_7.3-65           
[43] cli_3.6.6              mvtnorm_1.4-2          rmarkdown_2.31        
[46] rms_8.1-1              generics_0.1.4         otel_0.2.0            
[49] rstudioapi_0.19.0      future.apply_1.20.2    RcppArmadillo_15.4.0-1
[52] stringr_1.6.0          splines_4.6.1          parallel_4.6.1        
[55] BiocManager_1.30.27    base64enc_0.1-6        vctrs_0.7.3           
[58] sandwich_3.1-2         glmnet_5.0             Matrix_1.7-5          
[61] SparseM_1.84-2         callr_3.8.0            Formula_1.2-5         
[64] htmlTable_2.5.0        listenv_1.0.0          foreach_1.5.2         
[67] cmprsk_2.2-12          glue_1.8.1             parallelly_1.48.0     
[70] codetools_0.2-20       ps_1.9.3               stringi_1.8.7         
[73] shape_1.4.6.1          gtable_0.3.6           tibble_3.3.1          
[76] pillar_1.11.1          htmltools_0.5.9        quantreg_6.1          
[79] lava_1.9.2             R6_2.6.1               evaluate_1.0.5        
[82] lattice_0.22-9         backports_1.5.1        renv_1.2.3            
[85] MatrixModels_0.5-4     Rcpp_1.1.2             gridExtra_2.3.1       
[88] nlme_3.1-169           prodlim_2026.03.11     checkmate_2.3.4       
[91] xfun_0.60              zoo_1.8-15             pkgconfig_2.0.3       
```

## 4. Cohorte TCGA-LUAD

| Brazo | n | Eventos | Eventos % | Edad mediana [IQR] | Mujeres n (%) | Seguim. mediano (d) |
|---|---|---|---|---|---|---|
| Total | 504 | 182 | 36.1 | 66 [59-73] | 272 (54.0) | 882 |
| Training | 351 | 117 | 33.3 | 65 [58-72] | 185 (52.7) | 845 |
| Test | 153 | 65 | 42.5 | 67 [60-74] | 87 (56.9) | 1013 |

Seguimiento mediano por Kaplan-Meier inverso. Endpoint: **OS** (TCGA-CDR curado).

**Distribución por estadio (AJCC agrupado):**

| Estadio | Test | Training | Total |
|---|---|---|---|
| I | 81 | 189 | 270 |
| II | 36 | 84 | 120 |
| III | 25 | 56 | 81 |
| IV | 8 | 17 | 25 |
| Sin estadio | 3 | 5 | 8 |

**Casos completos usados en el modelado** (tras excluir NA en edad/sexo/estadio/score): training n=340 (116 eventos), test n=146 (61 eventos).

## 5. QC y preprocesamiento

- **Filtro de genes** (aprendido SOLO en training: >= 10 conteos en >= 20% de las muestras): **23770 de 60660 genes conservados** (39.2%).
- **Transformación**: vst (DESeq2; tendencia de dispersión ajustada en training y aplicada a todo).
- **QC de muestras** (|MAD z| > 4 en log10 tamaño de librería): **0 marcadas como atípicas de 504**. Ninguna descartada.

**Efecto de lote — R² de cada PC frente a la variable de lote (PCA sobre VST):**

| batch | pc | r2 |
|---|---|---|
| tss | PC1 | 0.106 |
| tss | PC2 | 0.063 |
| tss | PC3 | 0.159 |
| tss | PC4 | 0.070 |
| tss | PC5 | 0.052 |
| plate | PC1 | 0.059 |
| plate | PC2 | 0.042 |
| plate | PC3 | 0.069 |
| plate | PC4 | 0.050 |
| plate | PC5 | 0.047 |

Varianza explicada: PC1=9.7%, PC2=7.4%, PC3=5.9%, PC4=5.5%, PC5=3.9%.

> Nota: el lote se **documenta, no se corrige**. Una corrección agresiva puede eliminar señal biológica confundida con el centro (TSS) o la placa. Ver `results/figures/qc_pca_batch.png`.

## 6. Asociación univariable del score con la supervivencia (training)

| Termino | HR | IC95% inf | IC95% sup | p | C-index univariable |
|---|---|---|---|---|---|
| prolif_score (por DE) | 1.385 | 1.152 | 1.665 | 0.0005 | 0.614 |

Test log-rank global del modelo univariable: p = 0.0005 (n=340, eventos=116).

## 7. Modelos de Cox: clínico vs clínico + score (training)

| Modelo | Termino | HR | IC95% inf | IC95% sup | p |
|---|---|---|---|---|---|
| Clínico | age_std | 1.067 | 0.882 | 1.289 | 0.5056 |
| Clínico | gendermale | 1.115 | 0.769 | 1.618 | 0.5651 |
| Clínico | stage_groupII | 2.614 | 1.652 | 4.135 | 4.07e-05 |
| Clínico | stage_groupIII | 3.377 | 2.096 | 5.442 | 5.72e-07 |
| Clínico | stage_groupIV | 3.354 | 1.655 | 6.798 | 0.0008 |
| Clínico+Score | age_std | 1.152 | 0.945 | 1.403 | 0.1621 |
| Clínico+Score | gendermale | 1.090 | 0.753 | 1.579 | 0.6481 |
| Clínico+Score | stage_groupII | 2.455 | 1.547 | 3.895 | 0.0001 |
| Clínico+Score | stage_groupIII | 3.179 | 1.969 | 5.133 | 2.22e-06 |
| Clínico+Score | stage_groupIV | 3.001 | 1.470 | 6.125 | 0.0025 |
| Clínico+Score | prolif_score | 1.346 | 1.108 | 1.634 | 0.0027 |

Referencia de estadio: **I**. `age_std` está estandarizada con media/DE del **training**.

## 8. Test de razón de verosimilitud (modelos anidados)

| Modelo | loglik | Chisq | Df | p |
|---|---|---|---|---|
| Clínico | -556.5554 |    NA | NA | NA |
| Clínico+Score | -552.0045 | 9.1017 | 1 | 0.0026 |


## 9. Proporcionalidad de riesgos (Schoenfeld, modelo clínico+score)

| Termino | chisq | df | p | Flag |
|---|---|---|---|---|
| age_std | 1.7273 | 1 | 0.1888 | ok |
| gender | 2.5824 | 1 | 0.1081 | ok |
| stage_group | 6.2592 | 3 | 0.0997 | limitrofe (p<0.10) |
| prolif_score | 1.1780 | 1 | 0.2778 | ok |
| GLOBAL | 12.4819 | 6 | 0.0520 | limitrofe (p<0.10) |


## 10. Validación interna

### 10.1 C-index corregido por optimismo (bootstrap B=200)

| Modelo | C aparente | C corregido | Optimismo |
|---|---|---|---|
| Clinico | 0.6786 | 0.6628 | 0.0158 |
| Clinico+Score | 0.7028 | 0.6873 | 0.0155 |

**Incremento del C corregido (score vs clínico): 0.0245**.

> ALCANCE: `rms::validate` devuelve el índice corregido puntual, **no un IC**, así que **este ΔC concreto no lleva IC**. El contraste formal de si el score añade discriminación o precisión se hace con AUC y Brier en §10.6 (interno) y §11.5 (externo), que **sí** traen IC95% y p-valor.

### 10.2 C-index en el test interno retenido (con IC)

| modelo | C | lower | upper | se |
|---|---|---|---|---|
| Clinico | 0.6338 | 0.5544 | 0.7132 | 0.0405 |
| Clinico+Score | 0.6677 | 0.5878 | 0.7475 | 0.0407 |


### 10.3 AUC tiempo-dependiente (test interno)

| model | times | AUC | se | lower | upper |
|---|---|---|---|---|---|
| Clinico | 365 | 0.6499 | 0.0685 | 0.5156 | 0.7842 |
| Clinico | 1095 | 0.6760 | 0.0615 | 0.5555 | 0.7965 |
| Clinico | 1825 | 0.6520 | 0.0780 | 0.4992 | 0.8048 |
| Clinico+Score | 365 | 0.6976 | 0.0633 | 0.5735 | 0.8217 |
| Clinico+Score | 1095 | 0.7068 | 0.0574 | 0.5944 | 0.8193 |
| Clinico+Score | 1825 | 0.6901 | 0.1042 | 0.4858 | 0.8943 |


### 10.4 Brier e IPA (test interno)

| model | times | IPA | Brier | se | lower | upper |
|---|---|---|---|---|---|---|
| Null model | 365 | 0.0000 | 0.1116 | 0.0210 | 0.0703 | 0.1528 |
| Null model | 1095 | 0.0000 | 0.2422 | 0.0092 | 0.2242 | 0.2602 |
| Null model | 1825 | 0.0000 | 0.1824 | 0.0319 | 0.1200 | 0.2449 |
| Clinico | 365 | 0.0389 | 0.1072 | 0.0201 | 0.0679 | 0.1465 |
| Clinico | 1095 | 0.1106 | 0.2154 | 0.0200 | 0.1763 | 0.2546 |
| Clinico | 1825 | -0.2128 | 0.2213 | 0.0216 | 0.1790 | 0.2635 |
| Clinico+Score | 365 | 0.0541 | 0.1055 | 0.0197 | 0.0669 | 0.1442 |
| Clinico+Score | 1095 | 0.1189 | 0.2134 | 0.0208 | 0.1726 | 0.2542 |
| Clinico+Score | 1825 | -0.2563 | 0.2292 | 0.0244 | 0.1813 | 0.2771 |


> IPA = 1 - Brier/Brier(nulo KM). **IPA negativo significa peor que el modelo nulo.** Ver los valores a 1825 días.

### 10.5 Calibración

Curva de calibración del modelo clínico+score en el test interno: `results/figures/calibration_test.png` (riskRegression, método local, horizonte 1825 días).

### 10.6 CONTRASTES FORMALES entre modelos anidados (test interno)

Diferencia **clínico+score − clínico**, con IC95% y p (riskRegression, mismo dato, modelos anidados). Δ AUC > 0 favorece al score; **Δ Brier < 0 favorece al score**.

**Δ AUC:**

| Horizonte (d) | Delta AUC (score - clinico) | IC95% inf | IC95% sup | p |
|---|---|---|---|---|
| 365 | 0.0477 | -0.0205 | 0.1159 | 0.1708 |
| 1095 | 0.0308 | -0.0438 | 0.1054 | 0.4179 |
| 1825 | 0.0380 | -0.1052 | 0.1812 | 0.6027 |

**Δ Brier:**

| Horizonte (d) | Delta Brier (score - clinico) | IC95% inf | IC95% sup | p |
|---|---|---|---|---|
| 365 | -0.0017 | -0.0057 | 0.0023 | 0.3988 |
| 1095 | -0.0020 | -0.0187 | 0.0147 | 0.8142 |
| 1825 | 0.0079 | -0.0154 | 0.0313 | 0.5050 |

> **Conclusión formal (interno): el score NO añade discriminación ni precisión de forma estadísticamente detectable.** Todos los IC de Δ AUC y Δ Brier cruzan el cero y todos los p > 0.05. Con n=146 y 61 eventos el test interno tiene poca potencia para este contraste.

## 11. Validación externa (GEO)

- **Dataset**: GSE68465 (plataforma GPL96), cohorte quirúrgica de adenocarcinoma de pulmón.
- **n con OS válido**: 442; **eventos**: 236; genes mapeados a símbolo: 13237.
- **n usado con el modelo congelado** (casos completos + niveles válidos): 439.
- Estadio derivado de TNM patológico (aprox. AJCC 7ª, M0) porque el GSE no trae grupo AJCC.

### 11.1 Discriminación en GEO

| Modelo | C | IC95% inf | IC95% sup |
|---|---|---|---|
| Score solo | 0.6102 | 0.5728 | 0.6477 |
| Clínico congelado | 0.6953 | 0.6591 | 0.7315 |
| Clínico+Score congelado | 0.7080 | 0.6747 | 0.7413 |

**AUC tiempo-dependiente en GEO:**

| Horizonte_dias | AUC modelo congelado | AUC score solo |
|---|---|---|
| 365 | 0.7824 | 0.6683 |
| 1095 | 0.7740 | 0.6437 |
| 1825 | 0.7386 | 0.6154 |


### 11.2 Brier/IPA del modelo congelado en GEO (calibración transportada)

| model | times | IPA | Brier | se | lower | upper |
|---|---|---|---|---|---|---|
| Null model | 365 | 0.0000 | 0.0975 | 0.0117 | 0.0746 | 0.1204 |
| Null model | 1095 | 0.0000 | 0.2181 | 0.0080 | 0.2023 | 0.2338 |
| Null model | 1825 | 0.0000 | 0.2477 | 0.0024 | 0.2431 | 0.2524 |
| Clinico | 365 | 0.0764 | 0.0901 | 0.0105 | 0.0695 | 0.1106 |
| Clinico | 1095 | 0.1744 | 0.1800 | 0.0084 | 0.1636 | 0.1964 |
| Clinico | 1825 | 0.1392 | 0.2133 | 0.0088 | 0.1961 | 0.2304 |
| Clinico+Score | 365 | 0.0211 | 0.0954 | 0.0127 | 0.0706 | 0.1203 |
| Clinico+Score | 1095 | 0.0433 | 0.2086 | 0.0140 | 0.1812 | 0.2360 |
| Clinico+Score | 1825 | 0.0202 | 0.2427 | 0.0131 | 0.2171 | 0.2683 |


### 11.3 Replicación de la asociación (re-ajuste DENTRO de GEO)

| term | HR | lower | upper | p |
|---|---|---|---|---|
| age_std | 1.3463 | 1.1734 | 1.5446 | 2.23e-05 |
| gendermale | 1.2587 | 0.9682 | 1.6364 | 0.0856 |
| stage_groupII | 2.2881 | 1.6777 | 3.1207 | 1.72e-07 |
| stage_groupIII | 4.4669 | 3.2365 | 6.1651 | 8.66e-20 |
| stage_groupIV |    NA |    NA |    NA | NA |
| prolif_score | 1.2402 | 1.0825 | 1.4208 | 0.0019 |

Ajustado por: age_std, gender, stage_group · n=439 · eventos=235.

### 11.4 ¿Se degradó el desempeño respecto a lo interno?

| Métrica | Interno (test TCGA) | Externo (GEO) | Δ |
|---|---|---|---|
| C-index modelo clínico+score | 0.6677 | 0.7080 | 0.0403 |
| HR del score (ajustado) | 1.346 | 1.240 | atenuación |

**IPA del modelo clínico+score, interno vs externo:**

| Horizonte (d) | IPA interno (test TCGA) | IPA externo (GEO) | Δ (ext - int) |
|---|---|---|---|
| 365 | 0.0541 | 0.0211 | -0.0330 |
| 1095 | 0.1189 | 0.0433 | -0.0755 |
| 1825 | -0.2563 | 0.0202 | 0.2765 |

**Respuesta explícita, por métrica:**

- **Discriminación: NO se degradó.** C-index del modelo congelado en GEO = 0.708 vs 0.668 en el test interno (Δ = +0.040). El AUC(t) externo (0.78/0.77/0.74) también supera al interno (0.70/0.71/0.69).
- **Precisión (IPA): la comparación interno→externo NO es una inferencia válida.** El IPA interno a 1825 d es negativo (-0.256) por **inestabilidad de muestra pequeña**: el test interno tiene n=146 con 61 eventos y muy pocos sujetos aún en riesgo a 5 años, de modo que el Brier del nulo KM y el del modelo se estiman con enorme varianza. **No debe leerse como que el modelo 'mejora' al salir a GEO.** Las dos cohortes difieren en composición, seguimiento y plataforma; sus IPA no son comparables entre sí.
- **La señal fiable sobre calibración es la comparación DENTRO de GEO** (mismo dato, mismos sujetos, modelos anidados): ver el hallazgo negativo de abajo y los contrastes formales (§11.5).
- **Tamaño del efecto: SÍ se atenuó.** HR del score 1.35 (TCGA) → 1.24 (GEO).

**HALLAZGO NEGATIVO — dentro de GEO, añadir el score EMPEORA la precisión transportada:**

| Horizonte (d) | IPA clínico | IPA clínico+score | Δ (score - clínico) |
|---|---|---|---|
| 365 | 0.0764 | 0.0211 | -0.0553 |
| 1095 | 0.1744 | 0.0433 | -0.1310 |
| 1825 | 0.1392 | 0.0202 | -0.1190 |

En los tres horizontes el modelo **clínico solo** tiene mayor IPA (menor Brier) que el clínico+score al transportarlo congelado a GEO. **El C-index apunta en dirección CONTRARIA**: 0.6953 (clínico) vs 0.7080 (clínico+score), es decir +0.0127 a favor del score.

> Esto es una **disociación discriminación / calibración**: el score mejora levemente el *ordenamiento* de riesgo en GEO (C-index) pero degrada la *precisión de las probabilidades* predichas (Brier/IPA). Ambas cosas son compatibles y ambas deben reportarse.

> Lectura: la **discriminación** del score replica (HR 1.24 ajustado dentro de GEO, tabla 11.3), pero el **transporte congelado del modelo con score está peor calibrado** que el clínico. Causa probable: el score se estandariza con la media/DE del *training de TCGA* y esa escala no coincide con la distribución en GEO (plataforma microarray vs RNA-seq), de modo que el predictor lineal se desplaza. Esto es un problema de **recalibración**, no de ausencia de señal — pero tal como está, **el modelo con score no debe usarse congelado en una cohorte nueva sin recalibrar**.

> Interpretación cauta: que el desempeño externo sea *mejor* no indica un modelo superior en GEO, sino diferencias de composición de la cohorte (GEO es quirúrgica, sin estadio IV, con seguimiento más largo: mediana 1431 d vs 882 d en TCGA). El test interno de TCGA (n=146, 61 eventos) es pequeño y sus IC son anchos. **La comparación interno-externo no es un contraste formal.**

### 11.5 CONTRASTES FORMALES entre modelos anidados (GEO)

Diferencia **clínico+score − clínico** dentro de GEO (mismos sujetos). **Δ Brier > 0 significa que el score EMPEORA la precisión.**

**Δ AUC:**

| Horizonte (d) | Delta AUC (score - clinico) | IC95% inf | IC95% sup | p |
|---|---|---|---|---|
| 365 | 0.0409 | 0.0024 | 0.0795 | 0.0376 |
| 1095 | 0.0189 | -0.0135 | 0.0513 | 0.2531 |
| 1825 | 0.0134 | -0.0212 | 0.0480 | 0.4481 |

**Δ Brier:**

| Horizonte (d) | Delta Brier (score - clinico) | IC95% inf | IC95% sup | p |
|---|---|---|---|---|
| 365 | 0.0054 | 0.0004 | 0.0104 | 0.0357 |
| 1095 | 0.0286 | 0.0113 | 0.0458 | 0.0012 |
| 1825 | 0.0295 | 0.0053 | 0.0536 | 0.0167 |

> **Conclusión formal (externo): al transportar el modelo congelado, el score EMPEORA la precisión de forma estadísticamente significativa** (Δ Brier > 0, p < 0.05 en los tres horizontes), mientras que la discriminación (Δ AUC) no cambia de forma detectable. Ver §11.6: esto se debe a la ESCALA del score, y es corregible.

### 11.6 EXPERIMENTO: ¿la mala calibración externa es por la ESCALA del score?

**Diseño.** Se aplica el MISMO modelo congelado de TCGA a GEO en dos versiones que difieren *solo* en cómo se estandariza `prolif_score` (edad, sexo y estadio quedan idénticos):

- **(a) escala TCGA** — z-score con media/DE del *training* de TCGA (lo reportado hasta ahora).
- **(b) escala GEO** — z-score con media/DE de la propia cohorte GEO (recalibración de escala).

| Horizonte (d) | IPA_clinico | Brier_clinico | IPA_a_escalaTCGA | Brier_a_escalaTCGA | IPA_b_escalaGEO | Brier_b_escalaGEO |
|---|---|---|---|---|---|---|
| 365 | 0.0764 | 0.0901 | 0.0211 | 0.0954 | 0.1019 | 0.0876 |
| 1095 | 0.1744 | 0.1800 | 0.0433 | 0.2086 | 0.2136 | 0.1715 |
| 1825 | 0.1392 | 0.2133 | 0.0202 | 0.2427 | 0.1726 | 0.2050 |

**Resultado.** Re-escalar el score dentro de GEO **restaura la calibración y la mejora por encima del modelo clínico**: el IPA pasa de 0.021/0.043/0.020 (escala TCGA) a 0.102/0.214/0.173 (escala GEO), frente a 0.076/0.174/0.139 del clínico solo.

**Qué recalibra exactamente el re-escalado (descomposición algebraica).** La lectura NO es "escala del score *en lugar de* hazard basal": el re-escalado recalibra **dos cosas a la vez**. Con `z_T=(r-muT)/sdT` y `z_G=(r-muG)/sdG` se cumple la identidad

```
  beta * z_T  =  beta*(sdG/sdT) * z_G  +  beta*(muG - muT)/sdT
                 \_____pendiente_____/     \____constante (todos los sujetos)____/
```

- El **componente de MEDIA** desplaza el predictor lineal de **todos** los sujetos en una constante → esto es **calibración-en-lo-grande**, es decir **recalibración del basal/intercepto**.
- El **componente de DE** ajusta la **pendiente** del término del score.

| Componente | Valor | Debería ser |
|---|---|---|
| Media del score (GEO, en escala TCGA) | -3.179 | ~0 |
| → desplazamiento constante del lp | -0.944 | 0 |
| → factor multiplicativo del riesgo | 0.389 | 1 |
| DE del score (GEO, en escala TCGA) | 1.219 | ~1 |
| → pendiente efectiva del score | 0.362 | = beta = 0.297 |
| beta original (por DE de TCGA) | 0.297 | — |

> **Conclusión del experimento (corregida).** La descalibración externa es un problema de **recalibración con AMBOS ingredientes**, y el dominante es precisamente el del basal: en escala TCGA el score de GEO tiene media z = -3.18 (en vez de ~0), lo que desplaza el predictor lineal de todos los sujetos en -0.944 y multiplica el riesgo predicho por 0.389 — una infra-estimación sistemática de ~61%. A la vez, la DE (1.22 en vez de ~1) infla la pendiente efectiva de 0.297 a 0.362 (+22%). Re-estandarizar el score en la cohorte destino corrige **los dos simultáneamente**.

> **Implicación práctica (intacta):** el score es transportable **en rango** —el orden de riesgo se conserva, por eso la discriminación aguanta— pero **su estandarización debe re-estimarse en cada cohorte**; no debe usarse congelado el z-score de TCGA.

## 12. Score data-driven (elastic-net) frente al score simple y al clínico

Elastic-net Cox ajustado **solo en training** (glmnet, alpha=0.5, lambda por CV de 10 folds, prefiltro de 5000 genes de máxima varianza).

- **Genes seleccionados: 26** · lambda.min = 0.1695
- En GEO **solo 15 de los 26 genes están disponibles** (la plataforma GPL96 no cubre el resto), así que el predictor externo del elastic-net está truncado.

**Comparación lado a lado (C-index):**

| Modelo | C train aparente | C train corregido | C test interno | C GEO congelado |
|---|---|---|---|---|
| Clínico (edad+sexo+estadio) | 0.6786 | 0.6628 | 0.6338 | 0.6953 |
| Clínico + score proliferación | 0.7028 | 0.6873 | 0.6677 | 0.7080 |
| Score proliferación solo | 0.6140 | n/d | 0.5762 | 0.6102 |
| Elastic-net (26 genes) solo | 0.7740 | 0.7761 | 0.5637 | 0.5986 |

IC95% del C en el test interno: elastic-net 0.478–0.650 · score simple 0.499–0.654. En GEO: elastic-net 0.558–0.639.

> **Interpretación explícita: el modelo sofisticado NO supera al score simple; lo empeora.**
> - En training el elastic-net parece excelente (C=0.774) pero en el test interno **colapsa a 0.564** (IC 0.478–0.650, que **incluye 0.5**, es decir no distinguible del azar).
> - Es **peor que el modelo clínico solo** (0.634) y peor que clínico+score (0.668) en el mismo test.
> - La brecha aparente→test es de 0.210 puntos de C: **sobreajuste de manual**.
> - **AVISO METODOLÓGICO IMPORTANTE**: el 'C corregido' del elastic-net (0.776) es *engañoso* — `rms::validate` re-ajusta el Cox sobre el predictor YA seleccionado, **sin repetir la selección de genes dentro de cada bootstrap**, así que no captura el optimismo de la selección. Por eso da ≈ al aparente mientras el test real cae 0.21. Una corrección honesta exigiría envolver TODA la selección en el bootstrap.

## 13. Análisis de sensibilidad

### 13.1 Cox estratificado por estadio

| term | HR | lower | upper | p |
|---|---|---|---|---|
| age_std | 1.1549 | 0.9452 | 1.4111 | 0.1588 |
| gendermale | 1.0948 | 0.7529 | 1.5920 | 0.6353 |
| prolif_score | 1.3425 | 1.1006 | 1.6375 | 0.0037 |

**Re-chequeo de proporcionalidad (Schoenfeld) bajo estratificación:**

| Termino | chisq | df | p |
|---|---|---|---|
| age_std | 1.4428 | 1 | 0.2297 |
| gender | 2.0586 | 1 | 0.1513 |
| prolif_score | 1.1669 | 1 | 0.2800 |
| GLOBAL | 5.1413 | 3 | 0.1617 |

C-index dentro de estratos: clínico 0.5307 vs clínico+score 0.5894 (incremento 0.0588).

> Al estratificar por estadio, el HR del score se mantiene prácticamente idéntico (1.343 vs 1.346 en el modelo no estratificado): **la asociación no depende de la forma funcional del estadio**. Además el global de Schoenfeld deja de ser limítrofe (p=0.1617 vs 0.052 sin estratificar). NOTA: estos C-index (~0.53–0.59) NO son comparables con los de §10 porque el predictor lineal de un modelo estratificado **excluye** el efecto del estadio; miden solo la discriminación DENTRO de cada estrato.

### 13.2 Imputación múltiple (mice) vs caso completo

| Análisis | HR | IC95% inf | IC95% sup | p |
|---|---|---|---|---|
| Caso completo (primario) | 1.346 | 1.108 | 1.634 | 0.0027 |
| Imputación múltiple (mice, m=10) | 1.325 | 1.091 | 1.609 | 0.0049 |

> **Sin divergencia**: el HR del score es prácticamente el mismo con imputación múltiple que con caso completo. Los faltantes (8 muestras sin estadio) **no son una fuente de fragilidad** para esta conclusión.

### 13.3 RMST por grupo de score (training)

| Grupo | records | n.max | n.start | events | rmean | se(rmean) | median | 0.95LCL | 0.95UCL |
|---|---|---|---|---|---|---|---|---|---|
| score_group=score_alto | 170 | 170 | 170 | 70 | 1152.1 | 59.1 | 1229 | 864 | 2617 |
| score_group=score_bajo | 170 | 170 | 170 | 46 | 1397.0 | 55.5 | 1798 | 1268 | NA |

> Métrica libre del supuesto de proporcionalidad. Media de supervivencia restringida a 1825 días: **1152 d (score alto) vs 1397 d (score bajo)**, diferencia ~245 días a favor del score bajo. Consistente en dirección con el Cox.

## 14. Secciones Resultados y Discusión del reporte (verbatim)

# Resultados

## Descripción de la cohorte y partición

> _[chunk de código omitido — las cifras están en las secciones 4–11]_

## Control de calidad

> _[chunk de código omitido — las cifras están en las secciones 4–11]_

> _[chunk de código omitido — las cifras están en las secciones 4–11]_

## Exploración de supervivencia

> _[chunk de código omitido — las cifras están en las secciones 4–11]_

> _[chunk de código omitido — las cifras están en las secciones 4–11]_

> _[chunk de código omitido — las cifras están en las secciones 4–11]_

> _[chunk de código omitido — las cifras están en las secciones 4–11]_

## Modelo incremental (clínico vs clínico + score)

> _[chunk de código omitido — las cifras están en las secciones 4–11]_

> _[chunk de código omitido — las cifras están en las secciones 4–11]_

> **Lectura crítica.** El HR del score debe interpretarse **tras** el ajuste por
> estadio. Un score que solo es significativo sin ajuste es, muy probablemente,
> un proxy del estadio.

### Supuestos del modelo

> _[chunk de código omitido — las cifras están en las secciones 4–11]_

## Validación interna

> _[chunk de código omitido — las cifras están en las secciones 4–11]_

> _[chunk de código omitido — las cifras están en las secciones 4–11]_

> _[chunk de código omitido — las cifras están en las secciones 4–11]_

### Contrastes formales entre modelos anidados (test interno)

> _[chunk de código omitido — las cifras están en las secciones 4–11]_

> **Conclusión formal (interno):** todos los IC de Δ AUC y Δ Brier cruzan el cero y todos los
> p > 0.05 — el score **no añade discriminación ni precisión de forma detectable** en el test
> interno. Con 153 muestras y 65 eventos, este contraste tiene poca potencia.

## Score data-driven (elastic-net) frente al score simple

> _[chunk de código omitido — las cifras están en las secciones 4–11]_

> **El modelo sofisticado NO supera al score simple.** El elastic-net parece excelente en
> training (C≈0.77) pero **colapsa en el test interno (C≈0.56, IC incluye 0.5)**, por debajo
> del modelo clínico (0.63) y de clínico+score (0.67). La brecha aparente→test (~0.21) es
> sobreajuste de manual.
>
> **Aviso metodológico:** el "C corregido" del elastic-net es engañoso — `rms::validate`
> re-ajusta el Cox sobre un predictor **ya seleccionado**, sin repetir la selección de genes
> dentro de cada bootstrap, así que no captura el optimismo de la selección. Una corrección
> honesta exigiría envolver toda la selección en el bootstrap.

## Análisis de sensibilidad

> _[chunk de código omitido — las cifras están en las secciones 4–11]_

> **Sin divergencias que indiquen fragilidad.** El HR del score es estable al estratificar por
> estadio (≈1.34 vs 1.35) y con imputación múltiple (≈1.33 vs 1.35), y el global de Schoenfeld
> deja de ser limítrofe al estratificar. Los C-index del modelo estratificado (~0.53–0.59) **no**
> son comparables con los de la validación interna: el predictor lineal estratificado excluye el
> efecto del estadio y mide solo discriminación dentro de cada estrato.

## Validación externa (GEO)

> _[chunk de código omitido — las cifras están en las secciones 4–11]_

### Contrastes formales dentro de GEO

> _[chunk de código omitido — las cifras están en las secciones 4–11]_

> **Conclusión formal (externo):** al transportar el modelo **congelado**, el score
> **empeora la precisión de forma estadísticamente significativa** (Δ Brier > 0, p < 0.05 en
> los tres horizontes), mientras que la discriminación no cambia de forma detectable.

### Experimento: ¿la descalibración externa es por la ESCALA del score?

Se aplica el mismo modelo congelado a GEO en dos versiones que difieren *solo* en cómo se
estandariza `prolif_score` (edad, sexo y estadio idénticos): **(a)** z-score con media/DE del
training de TCGA, **(b)** z-score con media/DE de la propia cohorte GEO.

> _[chunk de código omitido — las cifras están en las secciones 4–11]_

> _[chunk de código omitido — las cifras están en las secciones 4–11]_

Re-estandarizar el score dentro de GEO restaura la calibración y la lleva **por encima** del
modelo clínico (IPA 0.10/0.21/0.17 frente a 0.08/0.17/0.14 del clínico), partiendo de
0.02/0.04/0.02 con la escala de TCGA.

**Qué recalibra exactamente.** La lectura **no** es "escala del score *en lugar de* hazard
basal": el re-escalado recalibra **dos cosas a la vez**. Con $z_T=(r-\mu_T)/\sigma_T$ y
$z_G=(r-\mu_G)/\sigma_G$:

$$\beta z_T \;=\; \underbrace{\beta\frac{\sigma_G}{\sigma_T}}_{\text{pendiente}} z_G \;+\; \underbrace{\beta\frac{\mu_G-\mu_T}{\sigma_T}}_{\text{constante para todos los sujetos}}$$

El componente de **media** desplaza el predictor lineal de *todos* los sujetos en una constante
— eso es **calibración-en-lo-grande**, es decir recalibración del **basal/intercepto**. El de
**DE** ajusta la **pendiente** del término del score.

> _[chunk de código omitido — las cifras están en las secciones 4–11]_

> **Conclusión corregida.** La descalibración externa es un problema de **recalibración con
> ambos ingredientes**, y el dominante es el del **basal**: en escala TCGA el score de GEO tiene
> media z ≈ −3.2 (en vez de ~0), lo que desplaza el predictor lineal de todos los sujetos y
> multiplica el riesgo predicho por ≈0.39 (infra-estimación sistemática de ~61%). A la vez la DE
> (≈1.22 en vez de ~1) infla la pendiente efectiva un ~22%. Re-estandarizar el score en la
> cohorte destino corrige **los dos simultáneamente**.
>
> **Implicación práctica (intacta):** el score es transportable **en rango** —el orden de riesgo
> se conserva, por eso la discriminación aguanta— pero **su estandarización debe re-estimarse en
> cada cohorte**; no debe usarse congelado el z-score de TCGA.

> **Lectura crítica.** La discriminación es invariante a re-escalados monótonos del score (por
> eso se transfiere), pero la **calibración** no lo es. Además, la comparación *interno →
> externo* del IPA **no es una inferencia válida**: el IPA interno a 5 años es negativo por
> inestabilidad de muestra pequeña (153 muestras, 65 eventos, muy pocos en riesgo a ese
> horizonte), no porque el modelo mejore al salir a GEO. La señal fiable sobre calibración es
> la comparación **dentro** de GEO, entre modelos anidados sobre los mismos sujetos.

# Discusión — amenazas a la validez

- **Potencia: el estudio mide asociación, no incremento de discriminación.** Con
  504 muestras y 182 eventos (153/65 en el test interno), el diseño está
  potenciado para detectar la **asociación** del score con la supervivencia
  (HR 1.35, p=0.003), pero **no** para estimar con precisión el **incremento de
  discriminación**: ΔC ≈ 0.01–0.03 con intervalos que se solapan y contrastes de
  Δ AUC / Δ Brier no significativos. Cuantificar ese incremento con utilidad
  clínica exigiría **muchos más eventos o el agrupamiento de varias cohortes**;
  con estos datos, afirmar una mejora de discriminación sería sobre-interpretar.
- **Ajuste por estadio (confusión dominante).** El estadio explica gran parte
  del pronóstico en LUAD; sin ajustarlo, cualquier firma proliferativa tiende a
  ser su proxy. Todo el reporte se lee *tras* ajuste.
- **El score data-driven no mejora al simple.** El elastic-net (26 genes) sobreajusta
  gravemente (C 0.77 en training → 0.56 en test) y queda por debajo del modelo
  clínico. La complejidad adicional no se traduce en desempeño.
- **Sesgo de espectro.** TCGA es una cohorte mayormente quirúrgica; la
  generalización a estadios avanzados o a población tratada es limitada.
- **Confusión residual.** Tabaquismo, tratamiento adyuvante y comorbilidad están
  mal capturados o ausentes; parte del "efecto" del score podría reflejarlos.
- **Ancestría.** TCGA está sesgada a ascendencia europea; la transportabilidad a
  otras poblaciones (p. ej. GSE31210, japonesa) no está garantizada.
- **Batch / plataforma.** Efecto de lote (TSS/placa) en TCGA y cambio de
  plataforma (RNA-seq → microarray) en GEO; se documenta y se prefiere un score
  rank-based en vez de sobre-corregir.
- **Multiplicidad y overfitting.** El elastic-net y la exploración por firmas
  invitan al sobreajuste; se mitiga con selección solo-en-training, corrección
  por optimismo y un test retenido.
- **Faltantes.** Análisis primario por caso completo con sensibilidad `mice`;
  divergencias entre ambos se interpretan como fragilidad.
- **Fiabilidad del endpoint.** Se usa TCGA-CDR (OS/PFI curados); aun así, la
  duración y calidad del seguimiento condicionan las estimaciones a largo plazo.


## 15. Código fuente (inline)

**Mapeo de nombres** — los archivos solicitados no existen con esos nombres; se incluyen los equivalentes reales por contenido:

| Solicitado | Real en el repo | Contenido |
|---|---|---|
| `R/04_scores_features.R` | `R/05_score.R` | singscore (G2M/E2F), escalado, elastic-net |
| `R/05_cox_models.R` | `R/06_model.R` | modelos Cox, supuestos |
| `R/06_internal_valid.R` | `R/06_model.R` | validación interna (mismo archivo) |
| `R/07_external_valid.R` | `R/07_external_valid.R` | validación externa GEO |

Se incluyen además `R/04_eda.R` (el 04 real) y `R/08_extra_analyses.R` (contrastes formales, evaluación del elastic-net y experimento de descalibración) para completitud.

### R/04_eda.R

```r
# =============================================================================
# 04_eda.R — Analisis exploratorio (descriptivo, sin decisiones de modelado)
# -----------------------------------------------------------------------------
# Describe la cohorte, el patron de faltantes y la supervivencia cruda global y
# por covariables. Es puramente descriptivo: NO se usa para seleccionar nada
# (la seleccion vive en training/05_score). El seguimiento mediano se estima por
# Kaplan-Meier inverso (metodo correcto: no es la mediana de 'time').
# =============================================================================

#' Tabla resumen de cohorte (estilo "Tabla 1").
cohort_summary_table <- function(clinical, params) {
  n  <- nrow(clinical)
  ev <- sum(clinical$status == 1L, na.rm = TRUE)
  # seguimiento mediano por KM inverso (censura como "evento")
  rkm <- survival::survfit(
    survival::Surv(time, status == 0L) ~ 1, data = clinical)
  medfu <- tryCatch(summary(rkm)$table[["median"]], error = function(e) NA_real_)
  stage_tab <- table(clinical$stage_group, useNA = "no")

  data.frame(
    variable = c("N muestras", "Eventos (n, %)", "Edad (mediana [IQR])",
                 "Mujeres (n, %)", "Estadio I/II/III/IV",
                 "Estadio faltante (n)", "Seguimiento mediano (dias, KM inverso)",
                 "Endpoint"),
    valor = c(
      as.character(n),
      sprintf("%d (%.1f%%)", ev, 100 * ev / n),
      sprintf("%.0f [%.0f-%.0f]",
              stats::median(clinical$age, na.rm = TRUE),
              stats::quantile(clinical$age, 0.25, na.rm = TRUE),
              stats::quantile(clinical$age, 0.75, na.rm = TRUE)),
      sprintf("%d (%.1f%%)", sum(clinical$gender == "female", na.rm = TRUE),
              100 * mean(clinical$gender == "female", na.rm = TRUE)),
      paste(as.integer(stage_tab[c("I", "II", "III", "IV")]), collapse = " / "),
      as.character(sum(is.na(clinical$stage_group))),
      as.character(round(medfu)),
      unique(clinical$endpoint)[1] %||% params$acquisition$endpoint
    ),
    stringsAsFactors = FALSE
  )
}

#' Mapa de faltantes por variable (proporcion de NA).
missingness_map <- function(clinical, params) {
  vars <- intersect(c("age", "gender", "stage_group", "time", "status",
                      "tss", "plate"), names(clinical))
  miss <- vapply(clinical[vars], function(x) mean(is.na(x)), numeric(1))
  data.frame(variable = names(miss),
             pct_missing = round(100 * miss, 2),
             row.names = NULL)
}

#' Grafico de barras del patron de faltantes.
plot_missingness <- function(miss_tab) {
  ggplot2::ggplot(miss_tab,
                  ggplot2::aes(stats::reorder(variable, pct_missing), pct_missing)) +
    ggplot2::geom_col(fill = "grey40") +
    ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = "% faltante", title = "Faltantes por variable") +
    ggplot2::theme_minimal()
}

#' Curva KM global (devuelve un ggplot para save_figure()).
km_global <- function(clinical, params) {
  fit <- survival::survfit(survival::Surv(time, status) ~ 1, data = clinical)
  p <- survminer::ggsurvplot(
    fit, data = clinical, conf.int = TRUE, risk.table = FALSE,
    xlab = "Dias", ylab = "Supervivencia",
    title = sprintf("KM global (%s)", unique(clinical$endpoint)[1]))
  p$plot
}

#' Curva KM estratificada por una variable categorica (estadio/sexo/edad-grupo).
km_by <- function(clinical, var, params) {
  df <- clinical
  if (var == "age_group") {
    df$age_group <- factor(ifelse(df$age >= 65, ">=65", "<65"),
                           levels = c("<65", ">=65"))
  }
  df <- df[!is.na(df[[var]]), , drop = FALSE]
  # Se incrusta el OBJETO formula en fit$call$formula: survminer con pval=TRUE
  # RE-EVALÚA fit$call$formula; si es as.formula(sprintf("... %s", var)), 'var'
  # se resuelve a la función base var() en su entorno -> error "coerce closure to
  # character". Con el objeto formula, re-evaluarlo devuelve la propia fórmula.
  f <- stats::as.formula(sprintf("Surv(time, status) ~ %s", var))
  fit <- survival::survfit(f, data = df)
  fit$call$formula <- f
  p <- survminer::ggsurvplot(
    fit, data = df, conf.int = FALSE, pval = TRUE, risk.table = FALSE,
    xlab = "Dias", ylab = "Supervivencia",
    title = sprintf("KM por %s", var), legend.title = var)
  p$plot
}
```

### R/05_score.R

```r
# =============================================================================
# 05_score.R — Score de proliferacion y features (SELECCION SOLO EN TRAINING)
# -----------------------------------------------------------------------------
# (a) Primario: score de proliferacion con singscore sobre Hallmark G2M/E2F.
#     singscore es rank-based POR MUESTRA: calcular el score en todas las
#     muestras a la vez NO es fuga (no ajusta nada entre muestras). Lo unico que
#     se aprende en training es la escala (z-score), que se aplica a test y GEO;
#     esto lo hace transferible entre plataformas sin re-entrenar.
# (b) Secundario: Cox elastic-net (glmnet) con lambda por CV, ajustado SOLO en
#     training. Prefiltro por varianza (training) para tratabilidad; el predictor
#     lineal se puede evaluar luego en test/GEO alineando genes.
# =============================================================================

#' Colapsa una matriz Ensembl x muestra a simbolo x muestra (maxima varianza).
expr_to_symbols <- function(expr_ens, gene_map) {
  gm <- gene_map[match(rownames(expr_ens), gene_map$ensembl), ]
  collapse_probes_to_symbols(expr_ens, gm$symbol)   # definida en 02_acquire.R
}

#' Descarga las firmas Hallmark solicitadas (simbolos) via msigdbr.
get_hallmark_sets <- function(params) {
  # msigdbr >= 10 usa 'collection' (antes 'category'); fallback por compatibilidad.
  md <- tryCatch(
    msigdbr::msigdbr(species = params$score$msigdb_species,
                     collection = params$score$msigdb_collection),
    error = function(e) msigdbr::msigdbr(species = params$score$msigdb_species,
                                         category = params$score$msigdb_collection))
  wanted <- params$score$gene_sets
  md <- md[md$gs_name %in% wanted, ]
  sets <- split(md$gene_symbol, md$gs_name)
  sets <- lapply(sets, unique)
  stopifnot(length(sets) >= 1)
  msg(sprintf("firmas Hallmark: %s",
              paste(sprintf("%s (%d genes)", names(sets), lengths(sets)),
                    collapse = ", ")))
  sets
}

#' Calcula el score de proliferacion con singscore. Devuelve el score combinado
#' (union de G2M+E2F como una unica firma 'up') y los scores por firma.
compute_singscore <- function(expr_sym, gene_sets, params) {
  rk <- singscore::rankGenes(expr_sym)
  prolif <- intersect(unique(unlist(gene_sets)), rownames(expr_sym))
  stopifnot(length(prolif) >= 10)
  sc <- singscore::simpleScore(rk, upSet = prolif)
  score <- stats::setNames(sc$TotalScore, colnames(expr_sym))

  per_set <- vapply(gene_sets, function(g) {
    g2 <- intersect(g, rownames(expr_sym))
    singscore::simpleScore(rk, upSet = g2)$TotalScore
  }, numeric(ncol(expr_sym)))
  rownames(per_set) <- colnames(expr_sym)

  list(score = score, per_set = as.data.frame(per_set), n_genes = length(prolif))
}

#' Aprende en training la escala del score y la aplica a todas las muestras.
#' scaling = 'zscore' (center/scale de training) o 'rank' (rango-normalizado).
scale_score <- function(score, train_ids, params) {
  if (identical(params$score$scaling, "rank")) {
    return(rank_normalize(score))   # global; robusto a plataforma
  }
  tr <- score[names(score) %in% train_ids]
  fit <- zscore_fit(tr)             # center/scale SOLO de training
  zscore_apply(score, fit)
}

#' Ensambla la tabla de scores para todas las muestras TCGA.
#' Devuelve data.frame: sample_id, prolif_raw, prolif_score (escalado), G2M, E2F.
assemble_scores <- function(vst_mat, gene_map, split, params) {
  expr_sym <- expr_to_symbols(vst_mat, gene_map)
  sets <- get_hallmark_sets(params)
  ss <- compute_singscore(expr_sym, sets, params)
  scaled <- scale_score(ss$score, split$train, params)

  df <- data.frame(
    sample_id   = names(ss$score),
    prolif_raw  = as.numeric(ss$score),
    prolif_score = as.numeric(scaled[names(ss$score)]),
    stringsAsFactors = FALSE
  )
  # anexar scores por firma con nombres cortos
  ps <- ss$per_set
  if ("HALLMARK_G2M_CHECKPOINT" %in% names(ps)) df$G2M <- ps[["HALLMARK_G2M_CHECKPOINT"]]
  if ("HALLMARK_E2F_TARGETS"   %in% names(ps)) df$E2F <- ps[["HALLMARK_E2F_TARGETS"]]
  df
}

# --- Score secundario: Cox elastic-net --------------------------------------

#' Ajusta Cox elastic-net en TRAINING. Prefiltra genes por varianza (training)
#' para tratabilidad. Devuelve el modelo, genes seleccionados y lambda.
fit_elasticnet_cox <- function(vst_mat, clinical, split, params) {
  common <- intersect(colnames(vst_mat), split$train)
  cl <- clinical[match(common, clinical$sample_id), ]
  ok <- is.finite(cl$time) & cl$time > 0 & cl$status %in% c(0L, 1L)
  common <- common[ok]; cl <- cl[ok, ]
  if (sum(cl$status == 1L) < params$score$elasticnet$min_events) {
    msg("elastic-net omitido: eventos insuficientes en training")
    return(NULL)
  }

  # prefiltro por varianza SOLO en training
  n_top <- params$score$elasticnet$top_var_genes %||% 5000
  Xtr <- vst_mat[, common, drop = FALSE]
  v <- matrixStats::rowVars(as.matrix(Xtr))
  top <- order(v, decreasing = TRUE)[seq_len(min(n_top, length(v)))]
  X <- t(Xtr[top, , drop = FALSE])            # muestras x genes

  set_global_seed(params)                     # CV reproducible
  y <- survival::Surv(cl$time, cl$status)
  cvfit <- glmnet::cv.glmnet(
    X, y, family = "cox",
    alpha = params$score$elasticnet$alpha,
    nfolds = params$score$elasticnet$nfolds)
  co <- as.matrix(stats::coef(cvfit, s = "lambda.min"))
  sel <- co[co[, 1] != 0, , drop = FALSE]
  msg(sprintf("elastic-net: %d genes seleccionados (lambda.min=%.4g)",
              nrow(sel), cvfit$lambda.min))
  list(cvfit = cvfit, coef = sel, genes = colnames(X),
       lambda_min = cvfit$lambda.min)
}

#' Predictor lineal del elastic-net en cualquier matriz genes x muestra,
#' alineando por nombre de gen (genes ausentes -> 0). Escalado en training.
predict_elasticnet <- function(fit, vst_mat, train_ids = NULL, params = NULL) {
  if (is.null(fit)) return(NULL)
  b <- fit$coef[, 1]
  g <- intersect(names(b), rownames(vst_mat))
  lp <- as.numeric(t(vst_mat[g, , drop = FALSE]) %*% b[g])
  names(lp) <- colnames(vst_mat)
  if (!is.null(train_ids) && !is.null(params)) {
    fit_sc <- zscore_fit(lp[names(lp) %in% train_ids])
    lp <- zscore_apply(lp, fit_sc)
  }
  lp
}
```

### R/06_model.R

```r
# =============================================================================
# 06_model.R — Modelado Cox incremental + validacion interna
# -----------------------------------------------------------------------------
# Compara el modelo CLINICO (edad, sexo, estadio) frente a CLINICO + score. El
# ajuste por estadio es obligatorio: sin el, una "firma" suele ser proxy del
# estadio. Se verifican supuestos (Schoenfeld, forma funcional, influencia) y se
# valida internamente con C-index corregido por optimismo (bootstrap), AUC
# tiempo-dependiente, calibracion y Brier/IPA. Analisis de sensibilidad para
# fallos de proporcionalidad y para faltantes (mice) y no linealidad (RMST/RSF).
# =============================================================================

# --- Ensamblado de datos de modelado ----------------------------------------

#' Une clinica + score, estandariza edad con parametros de TRAINING (sin fuga),
#' marca el brazo train/test y opcionalmente anexa el predictor elastic-net.
assemble_model_data <- function(clinical, scores, split, params, elastic_lp = NULL) {
  df <- merge(clinical,
              scores[, c("sample_id", "prolif_score", "prolif_raw")],
              by = "sample_id")
  tr <- df$sample_id %in% split$train
  age_fit <- zscore_fit(df$age[tr])              # center/scale de training
  df$age_std <- zscore_apply(df$age, age_fit)
  df$arm <- ifelse(df$sample_id %in% split$train, "train",
                   ifelse(df$sample_id %in% split$test, "test", NA_character_))
  df$stage_group <- droplevels(factor(df$stage_group,
                                      levels = c("I", "II", "III", "IV")))
  if (!is.null(elastic_lp)) df$elastic_lp <- elastic_lp[df$sample_id]
  df
}

#' Filtra a casos completos en las variables del modelo (analisis primario).
model_ready <- function(dat, params, extra = NULL) {
  vars <- unique(c("time", "status", params$model$clinical_vars,
                   "prolif_score", extra))
  vars <- intersect(vars, names(dat))
  dat[stats::complete.cases(dat[, vars, drop = FALSE]), , drop = FALSE]
}

.cox_formula <- function(rhs) {
  stats::as.formula(paste("survival::Surv(time, status) ~", rhs))
}

#' Modelo clinico (edad + sexo + estadio).
fit_cox_clinical <- function(dat, params) {
  rhs <- paste(params$model$clinical_vars, collapse = " + ")
  survival::coxph(.cox_formula(rhs), data = dat, x = TRUE, y = TRUE)
}

#' Modelo clinico + score de proliferacion.
fit_cox_full <- function(dat, params) {
  rhs <- paste(c(params$model$clinical_vars, "prolif_score"), collapse = " + ")
  survival::coxph(.cox_formula(rhs), data = dat, x = TRUE, y = TRUE)
}

#' Test de razon de verosimilitud clinico vs clinico+score (aporte incremental).
cox_lrt <- function(clin_fit, full_fit) {
  stats::anova(clin_fit, full_fit, test = "LRT")
}

#' Tabla ordenada de HR/IC/p de un coxph (para el reporte).
tidy_cox <- function(fit, ci = 0.95) {
  s <- summary(fit, conf.int = ci)
  cc <- s$conf.int
  co <- s$coefficients
  data.frame(
    term   = rownames(co),
    HR     = cc[, "exp(coef)"],
    lower  = cc[, 3],
    upper  = cc[, 4],
    p      = co[, "Pr(>|z|)"],
    row.names = NULL
  )
}

# --- Supuestos ---------------------------------------------------------------

#' Proporcionalidad de riesgos (Schoenfeld). Devuelve la tabla de cox.zph.
check_ph <- function(fit) {
  z <- survival::cox.zph(fit, transform = "km")
  as.data.frame(z$table)
}

#' Forma funcional del score: residuos de martingala de un modelo nulo frente al
#' score (una curva loess muy no lineal sugiere transformar/splines).
check_functional_form <- function(dat) {
  null <- survival::coxph(survival::Surv(time, status) ~ 1, data = dat)
  data.frame(prolif_score = dat$prolif_score,
             martingale   = as.numeric(stats::residuals(null, type = "martingale")))
}

#' Influencia por observacion (dfbeta) del modelo completo.
influence_dfbeta <- function(fit) {
  db <- stats::residuals(fit, type = "dfbeta")
  db <- as.data.frame(db)
  names(db) <- names(stats::coef(fit))
  db
}

# --- Discriminacion e integracion clinica -----------------------------------

#' Alinea 'newdata' a los niveles de factor vistos por el modelo (fit$xlevels) y
#' descarta filas con niveles no vistos (predict.coxph falla ante niveles nuevos;
#' esto ocurre con submuestras pequeñas o entre plataformas en validación externa).
align_to_model <- function(newdata, fit) {
  xl <- fit$xlevels
  if (is.null(xl)) return(newdata)
  keep <- rep(TRUE, nrow(newdata))
  for (v in intersect(names(xl), names(newdata))) {
    newdata[[v]] <- factor(as.character(newdata[[v]]), levels = xl[[v]])
    keep <- keep & !is.na(newdata[[v]])
  }
  newdata[keep, , drop = FALSE]
}

#' C-index de Harrell para un marcador (mayor marcador = mayor riesgo).
c_index <- function(time, status, marker, ci = 0.95) {
  cc <- survival::concordance(
    survival::Surv(time, status) ~ marker, reverse = TRUE)
  se <- sqrt(cc$var)
  z <- stats::qnorm(1 - (1 - ci) / 2)
  c(C = cc$concordance, lower = cc$concordance - z * se,
    upper = cc$concordance + z * se, se = se)
}

#' C-index corregido por optimismo via bootstrap (rms::validate -> Dxy).
#' C = Dxy/2 + 0.5. Reporta aparente y corregido.
optimism_corrected_c <- function(dat_train, rhs, params) {
  set_global_seed(params)
  dd <- rms::datadist(dat_train); on.exit(options(datadist = NULL), add = TRUE)
  options(datadist = dd)
  # rms::cph parsea mejor 'Surv' sin espacio de nombres (survival va adjunto en
  # los workers de targets vía tar_option_set(packages=...)).
  f <- stats::as.formula(paste("Surv(time, status) ~", rhs))
  fit <- rms::cph(f, data = dat_train, x = TRUE, y = TRUE, surv = TRUE)
  v <- rms::validate(fit, method = "boot", B = params$validation$bootstrap_B)
  dxy_app  <- v["Dxy", "index.orig"]
  dxy_corr <- v["Dxy", "index.corrected"]
  c(c_apparent = dxy_app / 2 + 0.5,
    c_optimism_corrected = dxy_corr / 2 + 0.5,
    optimism = (dxy_app - dxy_corr) / 2)
}

#' Horizontes válidos: descarta los que caen fuera del seguimiento observado
#' (no se puede evaluar discriminación/Brier más allá del último tiempo). Evita
#' errores con cohortes/submuestras de seguimiento corto (p. ej. modo smoke).
.valid_horizons <- function(time, params) {
  h <- params$validation$horizons_days
  h <- h[h < max(time, na.rm = TRUE)]
  if (length(h) == 0) h <- stats::median(time, na.rm = TRUE)
  h
}

#' AUC tiempo-dependiente (timeROC) de un marcador a los horizontes válidos.
td_auc <- function(time, status, marker, params) {
  timeROC::timeROC(T = time, delta = status, marker = marker, cause = 1,
                   times = .valid_horizons(time, params), iid = TRUE)
}

#' Metricas de riesgo (Brier/IPA + AUC) con censura tratada correctamente
#' (riskRegression::Score). Compara modelos en 'dat' a los horizontes de params.
score_metrics <- function(models, dat, params) {
  # OJO: formula con 'Surv' SIN espacio de nombres — riskRegression no reconoce
  # el tipo de respuesta si se usa survival::Surv (-> "Cannot assign response
  # type"). survival va adjunto en los workers de targets.
  riskRegression::Score(
    object   = models,                      # lista nombrada de coxph
    formula  = Surv(time, status) ~ 1,
    data     = dat,
    times    = .valid_horizons(dat$time, params),
    metrics  = c("AUC", "Brier"),
    summary  = "ipa",                        # Index of Prediction Accuracy
    conf.int = params$validation$ci_level,
    split.method = "none"
  )
}

#' Curva de calibracion a un horizonte dado (objeto de riskRegression::Score).
calibration_plot <- function(models, dat, params, horizon) {
  sc <- riskRegression::Score(
    object = models, formula = Surv(time, status) ~ 1,
    data = dat, times = horizon, plots = "calibration",
    conf.int = FALSE, split.method = "none")
  riskRegression::plotCalibration(sc, cens.method = "local", round = FALSE)
  invisible(sc)
}

# --- Analisis de sensibilidad -----------------------------------------------

#' Cox estratificado por estadio (si falla proporcionalidad del estadio).
fit_cox_stratified <- function(dat, params) {
  base <- setdiff(params$model$clinical_vars, "stage_group")
  rhs <- paste(c(base, "strata(stage_group)", "prolif_score"), collapse = " + ")
  survival::coxph(.cox_formula(rhs), data = dat, x = TRUE, y = TRUE)
}

#' RMST (media restringida) por grupo alto/bajo de score, hasta tau.
rmst_by_score <- function(dat, params, tau = NULL) {
  g <- factor(ifelse(dat$prolif_score >= stats::median(dat$prolif_score, na.rm = TRUE),
                     "score_alto", "score_bajo"))
  dat$score_group <- g
  if (is.null(tau)) tau <- max(params$validation$horizons_days)
  fit <- survival::survfit(survival::Surv(time, status) ~ score_group, data = dat)
  st <- summary(fit, rmean = tau)$table
  as.data.frame(st)
}

#' Cox completo con imputacion multiple (mice) como sensibilidad a faltantes.
cox_full_mice <- function(dat, params, m = 10) {
  vars <- intersect(c("time", "status", "age_std", "gender", "stage_group",
                      "prolif_score"), names(dat))
  sub <- dat[, vars, drop = FALSE]
  imp <- mice::mice(sub, m = m, printFlag = FALSE, seed = params$project$seed)
  rhs <- paste(c(params$model$clinical_vars, "prolif_score"), collapse = " + ")
  fits <- with(imp, survival::coxph(
    stats::as.formula(paste("survival::Surv(time, status) ~", rhs))))
  as.data.frame(summary(mice::pool(fits), conf.int = TRUE))
}

#' Random Survival Forest (sensibilidad no lineal). Opcional: requiere
#' randomForestSRC; si no esta instalado, se omite con aviso (no rompe el DAG).
fit_rsf <- function(dat, params) {
  if (!requireNamespace("randomForestSRC", quietly = TRUE)) {
    msg("RSF omitido: 'randomForestSRC' no instalado (analisis opcional).")
    return(NULL)
  }
  set_global_seed(params)
  vars <- intersect(c("age_std", "gender", "stage_group", "prolif_score"),
                    names(dat))
  f <- stats::as.formula(paste("Surv(time, status) ~", paste(vars, collapse = " + ")))
  randomForestSRC::rfsrc(f, data = dat, ntree = 500, importance = TRUE)
}
```

### R/07_external_valid.R

```r
# =============================================================================
# 07_external_valid.R — Validacion externa en cohorte GEO independiente
# -----------------------------------------------------------------------------
# Dos preguntas complementarias, honestas ante diferencias de plataforma:
#  (1) TRANSPORTE (modelo congelado): aplicar el Cox entrenado en TCGA a GEO y
#      recalcular discriminacion y calibracion. Se ESPERA degradacion por
#      plataforma/poblacion/estadificacion; la calibracion suele romperse mas
#      que la discriminacion (se discute recalibracion del hazard basal).
#  (2) REPLICACION de la asociacion: re-ajustar el Cox DENTRO de GEO con las
#      covariables disponibles y ver si el score sigue asociado tras ajuste.
# La discriminacion (C-index, AUC) es invariante a transformaciones monotonas
# del score, asi que el escalado no la afecta: por eso el score se transfiere.
# =============================================================================

# --- Escalados aprendidos en TRAINING de TCGA (para transporte congelado) ----

#' Escala z del score (parametros de TRAINING de TCGA) para transferir a GEO.
tcga_score_scaling <- function(scores, split) {
  zscore_fit(scores$prolif_raw[scores$sample_id %in% split$train])
}

#' Escala z de la edad (parametros de TRAINING de TCGA).
tcga_age_scaling <- function(clinical, split) {
  zscore_fit(clinical$age[clinical$sample_id %in% split$train])
}

# --- Score en GEO ------------------------------------------------------------

#' Score de proliferacion en GEO con las MISMAS firmas Hallmark y singscore.
#' Devuelve el score crudo (la escala se aplica en el ensamblado segun el uso).
geo_scores <- function(geo_cohort, params) {
  sets <- get_hallmark_sets(params)
  ss <- compute_singscore(geo_cohort$expr, sets, params)
  data.frame(sample_id = names(ss$score),
             prolif_raw = as.numeric(ss$score),
             stringsAsFactors = FALSE)
}

#' Datos GEO para RE-AJUSTE dentro de GEO (escala aprendida en GEO).
geo_refit_data <- function(geo_cohort, gscores, params) {
  df <- merge(geo_cohort$pheno, gscores, by = "sample_id")
  df$prolif_score <- zscore_apply(df$prolif_raw, zscore_fit(df$prolif_raw))
  df$age_std <- zscore_apply(df$age, zscore_fit(df$age))
  df$gender <- factor(as.character(df$gender), levels = c("female", "male"))
  df$stage_group <- factor(as.character(df$stage_group),
                           levels = c("I", "II", "III", "IV"))
  df
}

#' Datos GEO para TRANSPORTE congelado (escala aprendida en TCGA training).
geo_frozen_data <- function(geo_cohort, gscores, score_scaling, age_scaling, params) {
  df <- merge(geo_cohort$pheno, gscores, by = "sample_id")
  df$prolif_score <- zscore_apply(df$prolif_raw, score_scaling)
  df$age_std <- zscore_apply(df$age, age_scaling)
  df$gender <- factor(as.character(df$gender), levels = c("female", "male"))
  df$stage_group <- factor(as.character(df$stage_group),
                           levels = c("I", "II", "III", "IV"))
  df
}

# --- (1) Transporte: discriminacion del modelo congelado --------------------

#' Discriminacion en GEO: score solo y modelo TCGA congelado (C-index + AUC).
external_discrimination <- function(geo_frozen, tcga_full_fit, params) {
  ci <- params$validation$ci_level
  res <- list()
  res$score_C   <- c_index(geo_frozen$time, geo_frozen$status,
                           geo_frozen$prolif_raw, ci)
  res$score_auc <- td_auc(geo_frozen$time, geo_frozen$status,
                          geo_frozen$prolif_raw, params)

  # modelo congelado: predictor lineal sobre casos completos de sus terminos,
  # con niveles de factor alineados a los del modelo (descarta no vistos).
  terms <- all.vars(stats::formula(tcga_full_fit))
  terms <- intersect(terms, names(geo_frozen))
  cc <- geo_frozen[stats::complete.cases(geo_frozen[, terms, drop = FALSE]), ]
  cc <- align_to_model(cc, tcga_full_fit)
  lp <- as.numeric(stats::predict(tcga_full_fit, newdata = cc, type = "lp"))
  res$model_C   <- c_index(cc$time, cc$status, lp, ci)
  res$model_auc <- td_auc(cc$time, cc$status, lp, params)
  res$n_used <- nrow(cc)
  res
}

# --- (1) Transporte: calibracion y Brier/IPA en GEO -------------------------

#' Calibracion/Brier/IPA del modelo congelado en GEO (riskRegression::Score).
external_calibration <- function(tcga_clin_fit, tcga_full_fit, geo_frozen, params) {
  terms <- unique(c(all.vars(stats::formula(tcga_full_fit)),
                    all.vars(stats::formula(tcga_clin_fit))))
  terms <- intersect(terms, names(geo_frozen))
  cc <- geo_frozen[stats::complete.cases(geo_frozen[, terms, drop = FALSE]), ]
  cc <- align_to_model(cc, tcga_full_fit)   # full es superconjunto de los factores
  models <- list("Clinico" = tcga_clin_fit, "Clinico+Score" = tcga_full_fit)
  score_metrics(models, cc, params)
}

# --- (2) Replicacion: re-ajuste dentro de GEO -------------------------------

#' Re-ajusta el Cox dentro de GEO con las covariables disponibles (descarta las
#' ausentes o de un solo nivel) y reporta el HR del score tras ajuste.
external_refit_assoc <- function(geo_refit, params) {
  avail <- intersect(params$model$clinical_vars, names(geo_refit))
  avail <- avail[vapply(avail, function(v) {
    x <- geo_refit[[v]]; length(unique(x[!is.na(x)])) > 1
  }, logical(1))]
  keep <- stats::complete.cases(
    geo_refit[, c("time", "status", avail, "prolif_score"), drop = FALSE])
  dat <- geo_refit[keep, , drop = FALSE]
  rhs <- paste(c(avail, "prolif_score"), collapse = " + ")
  fit <- survival::coxph(.cox_formula(rhs), data = dat, x = TRUE, y = TRUE)
  list(table = tidy_cox(fit, params$validation$ci_level),
       adjusted_for = avail, n = nrow(dat),
       events = sum(dat$status == 1L))
}
```

### R/08_extra_analyses.R

```r
# =============================================================================
# 08_extra_analyses.R — Contrastes formales, elastic-net evaluado y experimento
#                       de descalibración externa
# -----------------------------------------------------------------------------
# Añade lo que estaba codificado pero sin reportar:
#  (1) contrastes formales AUC/Brier entre modelos anidados (riskRegression),
#  (2) desempeño del score data-driven (elastic-net) lado a lado con el simple,
#  (3) resumen del Cox estratificado por estadio,
#  (4) experimento que testea si la mala calibración externa se debe a la ESCALA
#      del score o al hazard basal transportado.
# =============================================================================

# --- (1) Contrastes formales -------------------------------------------------

#' Extrae los contrastes (modelo vs referencia) de un objeto riskRegression::Score.
#' Devuelve delta, IC95% y p por horizonte. NULL si el objeto no los trae.
contrast_table <- function(score_obj, metric = c("AUC", "Brier"),
                           model = "Clinico+Score", reference = "Clinico") {
  metric <- match.arg(metric)
  x <- score_obj[[metric]]$contrasts
  if (is.null(x)) return(NULL)
  d <- as.data.frame(x)
  d <- d[as.character(d$model) == model &
         as.character(d$reference) == reference, , drop = FALSE]
  dcol <- paste0("delta.", metric)
  if (!dcol %in% names(d)) dcol <- grep("^delta", names(d), value = TRUE)[1]
  d <- d[, c("times", dcol, "lower", "upper", "p"), drop = FALSE]
  names(d) <- c("Horizonte (d)", paste0("Delta ", metric, " (score - clinico)"),
                "IC95% inf", "IC95% sup", "p")
  rownames(d) <- NULL
  d
}

# --- (2) Elastic-net: predictor lineal y desempeño ---------------------------

#' Predictor lineal del elastic-net sobre una matriz con IDs Ensembl (TCGA).
elastic_lp_tcga <- function(elastic, vst) {
  b <- elastic$coef[, 1]
  g <- intersect(names(b), rownames(vst))
  stats::setNames(as.numeric(t(vst[g, , drop = FALSE]) %*% b[g]), colnames(vst))
}

#' Predictor lineal del elastic-net en GEO: mapea Ensembl -> símbolo. Los genes
#' del modelo ausentes en la plataforma GEO se pierden (se reporta cuántos).
elastic_lp_geo <- function(elastic, geo_expr, gene_map) {
  b <- elastic$coef[, 1]
  sym <- gene_map$symbol[match(names(b), gene_map$ensembl)]
  ok <- !is.na(sym) & sym %in% rownames(geo_expr)
  b2 <- tapply(b[ok], sym[ok], sum)          # colapsa símbolos duplicados
  g <- names(b2)
  list(lp = stats::setNames(as.numeric(t(geo_expr[g, , drop = FALSE]) %*% b2[g]),
                            colnames(geo_expr)),
       n_used = length(g), n_total = length(b))
}

#' Desempeño del elastic-net: C aparente y corregido por optimismo en training,
#' C en el test interno y C congelado en GEO.
#' ADVERTENCIA: la corrección por optimismo aquí SOLO corrige el reajuste de Cox
#' sobre un predictor YA seleccionado; NO reejecuta la selección de genes dentro
#' de cada bootstrap, así que SUBESTIMA el optimismo real del score data-driven.
elastic_performance <- function(elastic, vst, gene_map, moddata, geo_cohort, params) {
  if (is.null(elastic)) return(NULL)
  ci <- params$validation$ci_level
  d <- moddata
  d$elp <- elastic_lp_tcga(elastic, vst)[d$sample_id]
  tr <- d[d$arm == "train" & is.finite(d$elp) & is.finite(d$time), ]
  te <- d[d$arm == "test"  & is.finite(d$elp) & is.finite(d$time), ]

  set_global_seed(params)
  dd <- rms::datadist(tr); on.exit(options(datadist = NULL), add = TRUE)
  options(datadist = dd)
  fit <- rms::cph(Surv(time, status) ~ elp, data = tr, x = TRUE, y = TRUE, surv = TRUE)
  v <- rms::validate(fit, method = "boot", B = params$validation$bootstrap_B)

  g <- elastic_lp_geo(elastic, geo_cohort$expr, gene_map)
  ph <- geo_cohort$pheno
  ph$elp <- g$lp[ph$sample_id]
  ph <- ph[is.finite(ph$elp), ]

  list(n_genes = nrow(elastic$coef), lambda_min = elastic$lambda_min,
       c_apparent = v["Dxy", "index.orig"] / 2 + 0.5,
       c_corrected = v["Dxy", "index.corrected"] / 2 + 0.5,
       c_test = c_index(te$time, te$status, te$elp, ci),
       c_geo  = c_index(ph$time, ph$status, ph$elp, ci),
       geo_genes_used = g$n_used, geo_genes_total = g$n_total,
       n_train = nrow(tr), n_test = nrow(te), n_geo = nrow(ph))
}

# --- (3) Cox estratificado por estadio ---------------------------------------

#' Ajusta clínico y clínico+score ESTRATIFICADOS por estadio y devuelve HR del
#' score, re-chequeo de proporcionalidad y C-index de ambos (vía predictor lineal;
#' survival::concordance() falla sobre estos ajustes por el alcance de `data`).
strat_summary <- function(dat, params) {
  base <- setdiff(params$model$clinical_vars, "stage_group")
  f <- function(extra = NULL) stats::as.formula(paste(
    "Surv(time, status) ~", paste(c(base, "strata(stage_group)", extra), collapse = " + ")))
  fit_clin <- survival::coxph(f(), data = dat, x = TRUE, y = TRUE)
  fit_full <- survival::coxph(f("prolif_score"), data = dat, x = TRUE, y = TRUE)
  ci <- params$validation$ci_level
  lp_c <- as.numeric(stats::predict(fit_clin, newdata = dat, type = "lp"))
  lp_f <- as.numeric(stats::predict(fit_full, newdata = dat, type = "lp"))
  zz <- as.data.frame(survival::cox.zph(fit_full, transform = "km")$table)
  list(coef = tidy_cox(fit_full, ci),
       ph = cbind(Termino = rownames(zz), zz),
       c_clin = c_index(dat$time, dat$status, lp_c, ci),
       c_full = c_index(dat$time, dat$status, lp_f, ci))
}

# --- (4) Experimento de descalibración externa -------------------------------

#' Variante de geo_frozen con el score RE-ESTANDARIZADO con media/DE de la propia
#' cohorte GEO. Todo lo demás (edad, sexo, estadio) queda IDÉNTICO, para aislar
#' la hipótesis de que la mala calibración externa se debe a la escala del score.
geo_rescaled_data <- function(geo_frozen) {
  d <- geo_frozen
  d$prolif_score <- zscore_apply(d$prolif_raw, zscore_fit(d$prolif_raw))
  d
}

#' Descompone qué recalibra exactamente el re-escalado del score.
#' Con z_T = (r-muT)/sdT y z_G = (r-muG)/sdG se cumple:
#'     beta * z_T  =  beta*(sdG/sdT) * z_G  +  beta*(muG-muT)/sdT
#' es decir, el re-escalado toca DOS cosas a la vez:
#'  (i) COMPONENTE DE MEDIA -> desplaza el predictor lineal de TODOS los sujetos
#'      en una constante: es calibración-en-lo-grande, o sea recalibración del
#'      hazard basal / intercepto.
#'  (ii) COMPONENTE DE DE   -> cambia la PENDIENTE efectiva del término del score.
#' Por tanto NO es "escala en lugar de basal": es recalibración con ambos
#' ingredientes, y re-estandarizar en la cohorte destino los corrige a la vez.
rescale_decomposition <- function(sc_scaling, geo_frozen, cox_full) {
  b   <- stats::coef(cox_full)[["prolif_score"]]
  muT <- sc_scaling$center; sdT <- sc_scaling$scale
  muG <- mean(geo_frozen$prolif_raw, na.rm = TRUE)
  sdG <- stats::sd(geo_frozen$prolif_raw, na.rm = TRUE)
  off <- b * (muG - muT) / sdT
  list(beta = b, mu_tcga = muT, sd_tcga = sdT, mu_geo = muG, sd_geo = sdG,
       offset_lp = off, hr_factor = exp(off),
       slope_factor = sdG / sdT, beta_effective = b * (sdG / sdT),
       z_mean_geo = mean(geo_frozen$prolif_score, na.rm = TRUE),
       z_sd_geo   = stats::sd(geo_frozen$prolif_score, na.rm = TRUE))
}
```

