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
