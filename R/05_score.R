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
