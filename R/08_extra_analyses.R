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
