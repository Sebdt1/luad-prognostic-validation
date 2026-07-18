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
