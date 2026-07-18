# =============================================================================
# _targets.R — Orquestacion reproducible del pipeline (paquete targets)
# -----------------------------------------------------------------------------
# El DAG hace explicito el orden y las dependencias, y solo recomputa lo que
# cambia. Ejecutar:
#   targets::tar_make()            # corrida (usa run.smoke_test de params.yml)
#   targets::tar_visnetwork()      # ver el DAG
# El modo smoke se controla en config/params.yml (run.smoke_test: true/false).
# =============================================================================

library(targets)
library(tarchetypes)

# Funciones del proyecto (R/00_utils.R ... R/07_external_valid.R)
tar_source("R")

# Paquetes que targets carga en cada worker.
tar_option_set(
  packages = c(
    "here", "yaml",
    "SummarizedExperiment", "TCGAbiolinks", "GEOquery", "Biobase",
    "DESeq2", "edgeR", "limma", "matrixStats", "data.table",
    "msigdbr", "singscore",
    "glmnet", "survival", "survminer", "timeROC", "riskRegression", "rms",
    "mice", "readxl", "readr", "tidyr", "ggplot2", "patchwork"
  ),
  format = "rds"
)

list(

  # --- 0. Configuracion ------------------------------------------------------
  tar_target(params_file, here::here("config", "params.yml"), format = "file"),
  tar_target(params, load_params(params_file)),

  # --- 1-2. Adquisicion ------------------------------------------------------
  tar_target(se_tcga,      download_tcga_luad(params)),
  tar_target(counts_tcga,  tcga_counts(se_tcga)),
  tar_target(gene_map,     tcga_gene_map(se_tcga)),
  tar_target(cdr,          read_tcga_cdr(params)),
  tar_target(clinical_all, build_clinical(se_tcga, cdr, params)),

  # Submuestreo estratificado y determinista si smoke_test = TRUE.
  tar_target(clinical, {
    set_global_seed(params)
    if (isTRUE(params$run$smoke_test)) {
      ids <- maybe_subsample(clinical_all$sample_id, params,
                             strata = clinical_all$stage_group)
      clinical_all[clinical_all$sample_id %in% ids, , drop = FALSE]
    } else clinical_all
  }),

  tar_target(counts_sub,
    counts_tcga[, intersect(colnames(counts_tcga), clinical$sample_id),
                drop = FALSE]),

  # --- 1. Particion training/test (ANTES de todo analisis) -------------------
  tar_target(split, make_split(clinical, params)),
  tar_target(split_balance, split_balance_table(clinical, split, params)),

  # --- 3. QC / preprocesamiento ---------------------------------------------
  tar_target(keep_genes, run_gene_filter(counts_sub, split$train, params)),
  tar_target(vst,        vst_transform(counts_sub, keep_genes, split$train, params)),
  tar_target(sample_qc_tab, sample_qc(counts_sub, params)),
  tar_target(batch,      batch_pca(vst, clinical, params)),
  tar_target(fig_batch,  save_figure(plot_batch_pca(batch, "tss"),
                                     "qc_pca_batch", params), format = "file"),

  # --- 4. EDA ----------------------------------------------------------------
  tar_target(cohort_tab, cohort_summary_table(clinical, params)),
  tar_target(miss_tab,   missingness_map(clinical, params)),
  tar_target(fig_miss,   save_figure(plot_missingness(miss_tab),
                                     "eda_missingness", params), format = "file"),
  tar_target(fig_km,     save_figure(km_global(clinical, params),
                                     "eda_km_global", params), format = "file"),
  tar_target(fig_km_stage, save_figure(km_by(clinical, "stage_group", params),
                                       "eda_km_stage", params), format = "file"),

  # --- 5. Score / features (solo training para la seleccion) -----------------
  tar_target(scores,   assemble_scores(vst, gene_map, split, params)),
  tar_target(elastic,  fit_elasticnet_cox(vst, clinical, split, params)),

  # --- 6. Modelado Cox incremental ------------------------------------------
  tar_target(moddata,  assemble_model_data(clinical, scores, split, params)),
  tar_target(train_dat, model_ready(moddata[moddata$arm == "train", ], params)),
  tar_target(test_dat,  model_ready(moddata[moddata$arm == "test", ], params)),

  tar_target(cox_clin, fit_cox_clinical(train_dat, params)),
  tar_target(cox_full, fit_cox_full(train_dat, params)),
  tar_target(cox_clin_tab, tidy_cox(cox_clin, params$validation$ci_level)),
  tar_target(cox_full_tab, tidy_cox(cox_full, params$validation$ci_level)),
  tar_target(lrt,      as.data.frame(cox_lrt(cox_clin, cox_full))),

  # Supuestos
  tar_target(ph_tab,   check_ph(cox_full)),
  tar_target(ff_tab,   check_functional_form(train_dat)),
  tar_target(dfbeta_tab, influence_dfbeta(cox_full)),

  # --- 7. Validacion interna -------------------------------------------------
  tar_target(c_opt_clin, optimism_corrected_c(
    train_dat, paste(params$model$clinical_vars, collapse = " + "), params)),
  tar_target(c_opt_full, optimism_corrected_c(
    train_dat, paste(c(params$model$clinical_vars, "prolif_score"),
                     collapse = " + "), params)),
  tar_target(c_test, {
    td <- align_to_model(test_dat, cox_full)   # niveles de factor seguros
    lp_clin <- stats::predict(cox_clin, newdata = td, type = "lp")
    lp_full <- stats::predict(cox_full, newdata = td, type = "lp")
    rbind(
      data.frame(modelo = "Clinico",
                 as.list(c_index(td$time, td$status, lp_clin,
                                 params$validation$ci_level))),
      data.frame(modelo = "Clinico+Score",
                 as.list(c_index(td$time, td$status, lp_full,
                                 params$validation$ci_level)))
    )
  }),
  tar_target(metrics_test, score_metrics(
    list("Clinico" = cox_clin, "Clinico+Score" = cox_full), test_dat, params)),
  tar_target(fig_calib, {
    path <- file.path(here::here("results", "figures"), "calibration_test.png")
    grDevices::png(path, width = 1600, height = 1400, res = 200)
    calibration_plot(list("Clinico+Score" = cox_full), test_dat, params,
                     horizon = utils::tail(.valid_horizons(test_dat$time, params), 1))
    grDevices::dev.off()
    path
  }, format = "file"),

  # --- Sensibilidad ----------------------------------------------------------
  tar_target(cox_strat,  fit_cox_stratified(train_dat, params)),
  tar_target(rmst_tab,   rmst_by_score(train_dat, params)),
  tar_target(mice_tab,   cox_full_mice(moddata[moddata$arm == "train", ], params)),

  # --- 8. Validacion externa (GEO) ------------------------------------------
  tar_target(eset,        download_geo(params)),
  tar_target(geo_cohort,  build_geo_cohort(eset, params)),
  tar_target(gsc,         geo_scores(geo_cohort, params)),
  tar_target(sc_scaling,  tcga_score_scaling(scores, split)),
  tar_target(age_scaling, tcga_age_scaling(clinical, split)),
  tar_target(geo_frozen,  geo_frozen_data(geo_cohort, gsc, sc_scaling,
                                          age_scaling, params)),
  tar_target(geo_refit,   geo_refit_data(geo_cohort, gsc, params)),
  tar_target(ext_disc,    external_discrimination(geo_frozen, cox_full, params)),
  tar_target(ext_cal,     external_calibration(cox_clin, cox_full,
                                               geo_frozen, params)),
  tar_target(ext_assoc,   external_refit_assoc(geo_refit, params)),

  # --- 9. Contrastes formales, elastic-net evaluado y descalibración ---------
  tar_target(contr_int_auc,   contrast_table(metrics_test, "AUC")),
  tar_target(contr_int_brier, contrast_table(metrics_test, "Brier")),
  tar_target(contr_ext_auc,   contrast_table(ext_cal, "AUC")),
  tar_target(contr_ext_brier, contrast_table(ext_cal, "Brier")),
  tar_target(elastic_perf, elastic_performance(elastic, vst, gene_map, moddata,
                                               geo_cohort, params)),
  tar_target(strat_tab,    strat_summary(train_dat, params)),
  # Experimento: ¿la mala calibración externa es por la ESCALA del score?
  tar_target(geo_rescaled, geo_rescaled_data(geo_frozen)),
  tar_target(ext_cal_rescaled, external_calibration(cox_clin, cox_full,
                                                    geo_rescaled, params)),
  tar_target(rescale_decomp, rescale_decomposition(sc_scaling, geo_frozen, cox_full)),

  # --- Tablas exportadas para el reporte/deliverables ------------------------
  tar_target(t_cohort,  save_table(cohort_tab,  "cohort_summary", params), format = "file"),
  tar_target(t_balance, save_table(split_balance,"split_balance", params), format = "file"),
  tar_target(t_clin,    save_table(cox_clin_tab, "cox_clinical",  params), format = "file"),
  tar_target(t_full,    save_table(cox_full_tab, "cox_full",      params), format = "file"),
  tar_target(t_ctest,   save_table(c_test,       "cindex_test",   params), format = "file"),
  tar_target(t_ph,      save_table(cbind(term = rownames(ph_tab), ph_tab),
                                   "ph_schoenfeld", params), format = "file"),
  tar_target(t_extassoc, save_table(ext_assoc$table, "external_refit_assoc",
                                    params), format = "file"),

  # --- Reporte Quarto (autodetecta tar_read/tar_load usados en el .qmd) ------
  tar_quarto(report, path = "reports/report.qmd", quiet = FALSE)
)
