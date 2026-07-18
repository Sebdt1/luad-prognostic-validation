# =============================================================================
# 03_qc.R — Control de calidad y preprocesamiento
# -----------------------------------------------------------------------------
# - Filtro de genes de baja expresion: el umbral se APRENDE en training y se
#   aplica a todo (mismo conjunto de genes en train/test -> sin fuga).
# - Estabilizacion de varianza: VST (DESeq2) con la tendencia de dispersion
#   ajustada en training y aplicada al resto; alternativa logCPM (edgeR).
#   (La normalizacion no usa el desenlace: no es fuga selectiva, pero ajustar
#    en training aporta rigor extra y transferencia limpia.)
# - QC de muestras: atipicas por tamano de libreria (MADs sobre log10).
# - Deteccion de lote (TSS/placa) por PCA: se DOCUMENTA, no se sobre-corrige,
#   porque una correccion agresiva puede borrar senal biologica confundida.
# =============================================================================

#' Aprende en training el conjunto de genes a conservar (baja expresion fuera).
run_gene_filter <- function(counts, train_ids, params) {
  tr <- counts[, colnames(counts) %in% train_ids, drop = FALSE]
  keep <- rowMeans(tr >= params$qc$min_count) >= params$qc$min_prop_samples
  msg(sprintf("filtro de genes: %d/%d conservados (>=%d en >=%.0f%% del training)",
              sum(keep), length(keep), params$qc$min_count,
              100 * params$qc$min_prop_samples))
  keep
}

#' Transforma conteos a escala estabilizada. Devuelve matriz genes x muestras.
#' VST: tendencia de dispersion ajustada en training y transferida a todo.
vst_transform <- function(counts, keep, train_ids, params) {
  counts <- counts[keep, , drop = FALSE]

  if (identical(params$qc$transform, "logcpm")) {
    dge <- edgeR::DGEList(counts)
    dge <- edgeR::calcNormFactors(dge)          # TMM
    return(edgeR::cpm(dge, log = TRUE, prior.count = 1))
  }

  is_tr <- colnames(counts) %in% train_ids
  cd <- data.frame(
    row.names = colnames(counts),
    arm = ifelse(is_tr, "train", "test")
  )
  # 1) ajustar dispersion en training
  dds_tr <- DESeq2::DESeqDataSetFromMatrix(
    counts[, is_tr, drop = FALSE], colData = cd[is_tr, , drop = FALSE], design = ~1)
  dds_tr <- DESeq2::estimateSizeFactors(dds_tr)
  dds_tr <- DESeq2::estimateDispersions(dds_tr, fitType = "parametric")
  disp_fun <- DESeq2::dispersionFunction(dds_tr)

  # 2) aplicar la MISMA transformacion a todas las muestras
  dds_all <- DESeq2::DESeqDataSetFromMatrix(counts, colData = cd, design = ~1)
  dds_all <- DESeq2::estimateSizeFactors(dds_all)
  DESeq2::dispersionFunction(dds_all) <- disp_fun
  vsd <- DESeq2::varianceStabilizingTransformation(dds_all, blind = FALSE)
  SummarizedExperiment::assay(vsd)
}

#' QC de muestras: marca atipicas por tamano de libreria (nº de MADs, log10).
sample_qc <- function(counts, params) {
  libsize <- colSums(counts)
  lls <- log10(libsize + 1)
  med <- stats::median(lls)
  madv <- stats::mad(lls)
  z <- (lls - med) / (madv + 1e-9)
  out <- data.frame(
    sample_id     = names(libsize),
    library_size  = as.numeric(libsize),
    log10_libsize = lls,
    mad_z         = z,
    outlier       = abs(z) > params$qc$library_size_mad,
    stringsAsFactors = FALSE
  )
  msg(sprintf("QC muestras: %d atipicas por tamano de libreria (>%s MADs)",
              sum(out$outlier), params$qc$library_size_mad))
  out
}

#' PCA para detectar efecto de lote. Devuelve scores, varianza explicada y la
#' asociacion (R^2 de un lm de cada PC ~ variable de lote). NO corrige nada.
batch_pca <- function(vst_mat, clinical, params, n_top = 2000, n_pc = 5) {
  v <- matrixStats::rowVars(as.matrix(vst_mat))
  top <- order(v, decreasing = TRUE)[seq_len(min(n_top, length(v)))]
  pc <- stats::prcomp(t(vst_mat[top, , drop = FALSE]), scale. = FALSE)
  k <- min(n_pc, ncol(pc$x))
  scores <- as.data.frame(pc$x[, seq_len(k), drop = FALSE])
  scores$sample_id <- rownames(scores)

  keep_cols <- intersect(c("sample_id", "tss", "plate", "stage_group"),
                         names(clinical))
  merged <- merge(scores, clinical[, keep_cols, drop = FALSE], by = "sample_id")

  assoc <- do.call(rbind, lapply(params$qc$batch_vars, function(b) {
    if (!b %in% names(merged)) return(NULL)
    r2 <- vapply(paste0("PC", seq_len(k)), function(pcn) {
      tryCatch(summary(stats::lm(merged[[pcn]] ~ factor(merged[[b]])))$r.squared,
               error = function(e) NA_real_)
    }, numeric(1))
    data.frame(batch = b, pc = paste0("PC", seq_len(k)), r2 = as.numeric(r2))
  }))

  var_expl <- (pc$sdev^2 / sum(pc$sdev^2))[seq_len(k)]
  list(scores = merged, assoc = assoc, var_explained = var_expl)
}

#' Grafico PCA coloreado por una variable de lote (para el reporte).
plot_batch_pca <- function(batch, color_by = "tss") {
  df <- batch$scores
  ve <- batch$var_explained
  ggplot2::ggplot(df, ggplot2::aes(PC1, PC2, color = factor(.data[[color_by]]))) +
    ggplot2::geom_point(alpha = 0.7, size = 1.6) +
    ggplot2::labs(
      x = sprintf("PC1 (%.1f%%)", 100 * ve[1]),
      y = sprintf("PC2 (%.1f%%)", 100 * ve[2]),
      color = color_by,
      title = "PCA sobre expresion estabilizada (deteccion de lote)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")  # muchos niveles TSS: solo patron
}
