# =============================================================================
# dependencies.R — Manifiesto único de paquetes del proyecto
# -----------------------------------------------------------------------------
# Fuente de verdad de las dependencias. Lo consume:
#   - setup/bootstrap.R  -> las instala (CRAN + Bioconductor) y hace snapshot.
#   - renv::snapshot()   -> descubre paquetes por los library() de abajo.
# Mantener sincronizado con lo que realmente se usa en R/*.R.
# =============================================================================

# --- Clasificación por repositorio (para el instalador) ----------------------
# Lista alineada con el STACK del protocolo y con lo que el código realmente usa.
cran_pkgs <- c(
  "renv", "targets", "tarchetypes",
  "here", "yaml", "matrixStats", "data.table",
  "glmnet", "survival", "survminer", "timeROC", "riskRegression", "rms",
  "mice",
  "ggplot2", "patchwork", "tidyr", "readr", "readxl",
  "sessioninfo", "testthat", "knitr", "quarto"
)

bioc_pkgs <- c(
  "TCGAbiolinks", "SummarizedExperiment", "GEOquery", "Biobase",
  "DESeq2", "edgeR", "limma",
  "singscore", "GSVA", "msigdbr"
)

# --- Declaración para el descubrimiento de renv ------------------------------
# NB: envuelto en if(FALSE) para que el archivo sea 'source'-able sin cargar
# nada, pero renv::dependencies() sigue detectando estos library() por análisis
# estático del código fuente.
if (FALSE) {
  library(renv); library(targets); library(tarchetypes)
  library(here); library(yaml); library(matrixStats); library(data.table)
  library(glmnet); library(survival); library(survminer); library(timeROC)
  library(riskRegression); library(rms); library(mice)
  library(ggplot2); library(patchwork); library(tidyr); library(readr)
  library(readxl); library(sessioninfo); library(testthat); library(knitr)
  library(quarto)
  library(TCGAbiolinks); library(SummarizedExperiment); library(GEOquery)
  library(Biobase); library(DESeq2); library(edgeR); library(limma)
  library(singscore); library(GSVA); library(msigdbr)
}
