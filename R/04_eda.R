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
