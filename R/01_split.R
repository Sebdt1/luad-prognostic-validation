# =============================================================================
# 01_split.R — Partición training/test ANTES de cualquier análisis
# -----------------------------------------------------------------------------
# Requisito crítico anti-fuga: la partición se fija primero y TODA selección de
# features/score/umbral se hará solo en 'train'. El 'test' queda intacto como
# validación interna honesta. La partición:
#   - es determinista (semilla global),
#   - se estratifica por estadio (confusor pronóstico dominante) para que ambos
#     brazos tengan pronóstico basal comparable,
#   - conserva las muestras con estadio faltante como su propio estrato (no se
#     descartan aquí; el manejo de faltantes se decide en el modelado).
# La estratificación usa SOLO el estratificador prefijado, no el desenlace.
# =============================================================================

#' Construye la partición training/test estratificada por estadio.
#'
#' @param clinical data.frame con al menos: sample_id y la columna estratificadora
#'   indicada en params$split$stratify_by.
#' @param params lista de configuración (config/params.yml).
#' @return lista con vectores 'train' y 'test' de sample_id, más metadatos.
make_split <- function(clinical, params) {
  set_global_seed(params)                      # partición reproducible
  p <- as.numeric(params$split$p_train)
  strat_col <- params$split$stratify_by
  stopifnot(all(c("sample_id", strat_col) %in% names(clinical)))

  ids <- as.character(clinical$sample_id)
  strata <- as.character(clinical[[strat_col]])
  strata[is.na(strata) | strata == ""] <- "__NA__"   # estrato propio para NA

  train <- character(0)
  for (g in unique(strata)) {
    gi <- ids[strata == g]
    n_train <- floor(length(gi) * p)
    # al menos 1 por estrato en cada brazo si el estrato tiene >= 2
    if (length(gi) >= 2) n_train <- min(max(n_train, 1L), length(gi) - 1L)
    train <- c(train, sample(gi, n_train))
  }
  test <- setdiff(ids, train)

  res <- list(
    train = sort(train),
    test  = sort(test),
    p_train = p,
    stratify_by = strat_col,
    seed = as.integer(params$project$seed),
    n_train = length(train),
    n_test = length(test)
  )
  msg(sprintf("split: %d train / %d test (estrato: %s, p=%.2f)",
              res$n_train, res$n_test, strat_col, p))
  res
}

#' Tabla resumen del balance de estratos entre brazos (para el reporte/QC).
split_balance_table <- function(clinical, split, params) {
  strat_col <- params$split$stratify_by
  df <- clinical
  df$arm <- ifelse(df$sample_id %in% split$train, "train",
                   ifelse(df$sample_id %in% split$test, "test", NA))
  df <- df[!is.na(df$arm), ]
  tab <- as.data.frame(table(
    stratum = addNA(df[[strat_col]]),
    arm = df$arm
  ))
  tidyr::pivot_wider(tab, names_from = arm, values_from = Freq)
}

#' Restringe cualquier tabla/matriz a las muestras de un brazo dado.
#' Uso: subset_to(expr_mat, split$train) o subset_to(clinical, split$train).
subset_to <- function(x, ids) {
  if (is.data.frame(x)) {
    return(x[as.character(x$sample_id) %in% ids, , drop = FALSE])
  }
  # matriz genes x muestras -> columnas
  x[, colnames(x) %in% ids, drop = FALSE]
}
