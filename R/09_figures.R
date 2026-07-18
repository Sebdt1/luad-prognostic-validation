# =============================================================================
# 09_figures.R — Figuras del reporte
# -----------------------------------------------------------------------------
# Las figuras se GENERAN dentro de los chunks del .qmd (no se enlazan desde
# results/figures/): así Quarto las numera, les pone pie, permite referencias
# cruzadas y las embebe en el HTML autocontenido. Enlazar por ruta relativa
# fallaba en silencio porque el directorio de ejecución de Quarto no coincide
# con el del .qmd.
# =============================================================================

#' KM por grupo alto/bajo de score (mediana como corte; solo para visualización,
#' el modelado usa el score continuo).
plot_km_score <- function(dat, params) {
  d <- dat
  d$score_group <- factor(
    ifelse(d$prolif_score >= stats::median(d$prolif_score, na.rm = TRUE),
           "Score alto", "Score bajo"),
    levels = c("Score bajo", "Score alto"))
  km_by(d, "score_group", params)
}

#' Forest plot de HR con IC95% en escala logarítmica.
plot_forest <- function(tab, title = NULL) {
  d <- tab
  d$term <- factor(d$term, levels = rev(d$term))
  ggplot2::ggplot(d, ggplot2::aes(x = HR, y = term)) +
    ggplot2::geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = lower, xmax = upper), height = 0.18) +
    ggplot2::geom_point(size = 2.6) +
    ggplot2::scale_x_log10() +
    ggplot2::labs(x = "HR (IC95%, escala log)", y = NULL, title = title) +
    ggplot2::theme_minimal(base_size = 11)
}

#' AUC tiempo-dependiente por horizonte, con IC95%, para los modelos comparados.
plot_auc_time <- function(metrics) {
  d <- as.data.frame(metrics$AUC$score)
  d <- d[as.character(d$model) != "Null model", ]
  d$model <- as.character(d$model)
  pd <- ggplot2::position_dodge(width = 90)
  ggplot2::ggplot(d, ggplot2::aes(times, AUC, colour = model, group = model)) +
    ggplot2::geom_hline(yintercept = 0.5, linetype = 2, colour = "grey60") +
    ggplot2::geom_line(position = pd, linewidth = 0.7) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper),
                           width = 60, position = pd) +
    ggplot2::geom_point(size = 2.6, position = pd) +
    ggplot2::labs(x = "Horizonte (días)", y = "AUC(t) con IC95%", colour = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "top")
}

#' Recalibración en GEO: IPA del modelo clínico frente al clínico+score con la
#' escala congelada de TCGA y con la escala re-estimada en GEO.
plot_geo_recalib <- function(ext_cal, ext_cal_rescaled) {
  g <- function(sc, m, lab) {
    d <- as.data.frame(sc$Brier$score)
    d <- d[as.character(d$model) == m, c("times", "IPA")]
    d$serie <- lab; d
  }
  d <- rbind(g(ext_cal, "Clinico", "Clínico (referencia)"),
             g(ext_cal, "Clinico+Score", "Clínico+score — escala TCGA (congelada)"),
             g(ext_cal_rescaled, "Clinico+Score", "Clínico+score — escala re-estimada en GEO"))
  ggplot2::ggplot(d, ggplot2::aes(times, IPA, colour = serie, group = serie)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2.6) +
    ggplot2::labs(x = "Horizonte (días)", y = "IPA (mayor = mejor)", colour = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "top",
                   legend.direction = "vertical", legend.text = ggplot2::element_text(size = 8))
}
