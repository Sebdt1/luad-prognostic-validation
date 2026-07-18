# =============================================================================
# bootstrap.R — Instala el entorno reproducible y GENERA renv.lock
# -----------------------------------------------------------------------------
# Por qué un script y no un renv.lock escrito a mano: un lockfile válido debe
# contener versiones y hashes exactos capturados de un repositorio real. Un
# lockfile inventado no restaura de forma fiable. Aquí se instala el stack y se
# captura el estado real con renv::snapshot(), produciendo un lockfile honesto.
#
# Uso (desde la raíz del proyecto, una sola vez):
#   Rscript setup/bootstrap.R
# Después, cualquiera reproduce con:
#   renv::restore()
# =============================================================================

message(">> Bootstrap del entorno reproducible (renv + Bioconductor)")

# --- 1. renv -----------------------------------------------------------------
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv", repos = "https://cloud.r-project.org")
}

# Inicializa renv sin instalar aún (bare) si el proyecto no está bajo renv.
if (!file.exists("renv/activate.R")) {
  renv::init(bare = TRUE, restart = FALSE)
}

# --- 2. Repositorios: CRAN + Bioconductor ------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  renv::install("BiocManager")
}
# Alinea la versión de Bioconductor con la de R en uso.
bioc_version <- BiocManager::version()
message(sprintf(">> Bioconductor %s (R %s)", bioc_version, getRversion()))
options(repos = BiocManager::repositories())

# --- 3. Instalar dependencias del manifiesto ---------------------------------
source("setup/dependencies.R")

message(">> Instalando paquetes CRAN...")
renv::install(cran_pkgs)

message(">> Instalando paquetes Bioconductor...")
# renv enruta a Bioconductor con el prefijo 'bioc::'
renv::install(paste0("bioc::", bioc_pkgs))

# --- 4. Capturar el estado real en el lockfile -------------------------------
message(">> Generando renv.lock (snapshot del estado instalado)...")
renv::snapshot(prompt = FALSE)

message(">> Listo. Verifica con: renv::status()  y luego  targets::tar_make()")
