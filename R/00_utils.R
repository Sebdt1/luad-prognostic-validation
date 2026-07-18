# =============================================================================
# 00_utils.R — Utilidades transversales del pipeline
# -----------------------------------------------------------------------------
# Funciones puras y auxiliares (sin efectos secundarios ocultos) usadas por el
# resto de etapas. Se cargan por _targets.R vía tar_source(). Racional de cada
# helper en su cabecera. Los tests de tests/test_utils.R fijan su comportamiento.
# =============================================================================

# --- Infraestructura ---------------------------------------------------------

#' Operador de coalescencia de nulos/vacíos
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

#' Mensaje con marca de tiempo (trazabilidad en el log de la corrida)
msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")
  invisible(NULL)
}

#' Carga la configuración central (única fuente de verdad de parámetros)
load_params <- function(path = here::here("config", "params.yml")) {
  stopifnot(file.exists(path))
  yaml::read_yaml(path)
}

#' Fija la semilla global desde params (reproducibilidad; se reporta en el .qmd)
set_global_seed <- function(params) {
  seed <- as.integer(params$project$seed)
  set.seed(seed)
  # Semilla para paralelismo tipo L'Ecuyer si se usa
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  invisible(seed)
}

# --- Persistencia de artefactos ---------------------------------------------

#' Guarda una tabla en results/tables (CSV) y devuelve la ruta (para targets)
save_table <- function(x, name, params, dir = here::here("results", "tables")) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(dir, paste0(name, ".csv"))
  readr::write_csv(as.data.frame(x), path)
  msg("tabla escrita: ", path)
  path
}

#' Guarda una figura en results/figures y devuelve la ruta (para targets)
save_figure <- function(plot, name, params, width = 7, height = 5, dpi = 300,
                        dir = here::here("results", "figures")) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(dir, paste0(name, ".png"))
  ggplot2::ggsave(path, plot = plot, width = width, height = height,
                  dpi = dpi, bg = "white")
  msg("figura escrita: ", path)
  path
}

#' Caché idempotente de una expresión costosa (descargas). Si existe el RDS y
#' overwrite=FALSE, lo lee; si no, evalúa 'expr', lo guarda y lo devuelve.
#' (targets ya cachea entre targets; esto protege descargas dentro de un target.)
cache_rds <- function(expr, path, overwrite = FALSE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(path) && !overwrite) {
    msg("caché HIT: ", path)
    return(readRDS(path))
  }
  msg("caché MISS: computando -> ", path)
  obj <- force(expr)
  saveRDS(obj, path)
  obj
}

# --- Muestreo para 'smoke test' ---------------------------------------------

#' Submuestrea identificadores en modo smoke, de forma estratificada y
#' determinista (respeta la semilla ya fijada). Devuelve todos los ids si no
#' está activo el smoke o si hay menos ids que el objetivo.
maybe_subsample <- function(ids, params, strata = NULL) {
  if (!isTRUE(params$run$smoke_test)) return(ids)
  n_target <- as.integer(params$run$smoke_n)
  if (length(ids) <= n_target) return(ids)
  if (is.null(strata)) {
    return(sort(sample(ids, n_target)))
  }
  # Muestreo proporcional por estrato
  strata <- as.character(strata)
  keep <- unlist(lapply(split(ids, strata), function(g) {
    k <- max(1L, round(n_target * length(g) / length(ids)))
    sample(g, min(k, length(g)))
  }), use.names = FALSE)
  sort(keep)
}

# --- Armonización clínica ----------------------------------------------------

#' Colapsa el estadio patológico AJCC (cadenas TCGA-CDR) a I/II/III/IV.
#' Los valores no informativos ("[Not Available]", "[Discrepancy]", "", NA,
#' "Stage X", "I/II NOS") se devuelven como NA para tratarlos como faltantes.
#' Racional: el estadio fino (IA/IB...) fragmenta y no aporta señal pronóstica
#' estable con estas n; el grupo I-IV es el estándar clínico defendible.
collapse_stage <- function(stage_chr) {
  s <- toupper(trimws(as.character(stage_chr)))
  s <- gsub("^STAGE\\s+", "", s)
  out <- rep(NA_character_, length(s))
  out[grepl("^IV", s)]                 <- "IV"
  out[grepl("^III", s)]                <- "III"
  out[grepl("^II($|[AB])", s)]         <- "II"
  out[grepl("^I($|[AB])", s)]          <- "I"
  # 'IS', 'X', 'NOS', 'TIS' y no informativos quedan NA
  factor(out, levels = c("I", "II", "III", "IV"))
}

#' Deriva el grupo de estadio AJCC (I-IV) desde TNM patológico tipo "pN0pT1"
#' (formato de algunas cohortes GEO, p. ej. GSE68465). Aproximación AJCC 7ª ed.
#' asumiendo M0 (cohorte resecada). Documentar como aproximación en el reporte.
#' Reglas: N0 -> I (T1-2), II (T3), III (T4); N1 -> II (T1-2), III (T3-4);
#' N2+ -> III. NX/valores no informativos -> NA.
tnm_to_stage_group <- function(x) {
  s  <- toupper(as.character(x))
  Tn <- suppressWarnings(as.integer(sub(".*PT([0-9X]).*", "\\1", s)))
  Nn <- suppressWarnings(as.integer(sub(".*PN([0-9X]).*", "\\1", s)))
  out <- rep(NA_character_, length(s))
  ok <- !is.na(Tn) & !is.na(Nn)
  out[ok & Nn == 0 & Tn <= 2] <- "I"
  out[ok & Nn == 0 & Tn == 3] <- "II"
  out[ok & Nn == 0 & Tn >= 4] <- "III"
  out[ok & Nn == 1 & Tn <= 2] <- "II"
  out[ok & Nn == 1 & Tn >= 3] <- "III"
  out[ok & Nn >= 2]           <- "III"
  factor(out, levels = c("I", "II", "III", "IV"))
}

#' Estandariza un vector numérico guardando center/scale (para transferir la
#' MISMA escala aprendida en training a test y a la cohorte externa -> sin fuga).
zscore_fit <- function(x) {
  center <- mean(x, na.rm = TRUE)
  scale  <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(scale) || scale == 0) scale <- 1
  list(center = center, scale = scale)
}

#' Aplica una escala aprendida (de zscore_fit) a un vector nuevo.
zscore_apply <- function(x, fit) {
  (x - fit$center) / fit$scale
}

#' Rango-normaliza a [0,1] (alternativa robusta a plataforma para transferencia).
rank_normalize <- function(x) {
  r <- rank(x, na.last = "keep", ties.method = "average")
  (r - 1) / (sum(!is.na(x)) - 1)
}
