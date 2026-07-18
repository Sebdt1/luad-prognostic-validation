# =============================================================================
# 02_acquire.R — Adquisición de datos (idempotente y cacheada)
# -----------------------------------------------------------------------------
# Primario  : TCGA-LUAD RNA-seq (STAR counts) vía TCGAbiolinks (open-access).
# Endpoint  : recurso clínico curado TCGA-CDR (OS/PFI), NO days_to_death crudo.
# Externo   : cohorte GEO de LUAD con supervivencia vía GEOquery.
# Todas las descargas se cachean en data/raw/ y se saltan si ya existen.
# =============================================================================

# --- TCGA-LUAD RNA-seq -------------------------------------------------------

#' Descarga y prepara el SummarizedExperiment de STAR counts de TCGA-LUAD.
#' Cacheado en data/raw/tcga_luad_se.rds (idempotente: descarga UNA vez).
#' Se descarga la cohorte completa (el filtro `barcode=` de GDCquery es frágil a
#' escala); el submuestreo del modo smoke se aplica aguas abajo (target
#' `clinical`), así que el cómputo posterior sigue siendo rápido.
# -----------------------------------------------------------------------------
# NOTA IMPORTANTE (workaround): la descarga de RNA-seq de TCGAbiolinks 2.40.0
# está ROTA en este entorno por DOS bugs distintos:
#   - método "api": GDCdownload.aux() termina en message() (retorna NULL) y el
#     llamador hace `if (ret == 1)` -> error "argumento de longitud cero".
#   - método "client": descarga el zip de gdc-client pero no lo extrae/ejecuta.
# Solución robusta y reproducible: usamos GDCquery SOLO para el catálogo de
# archivos, descargamos con gdc-client (la herramienta oficial del GDC) mediante
# un manifiesto, y construimos el SummarizedExperiment leyendo los .tsv STAR
# directamente. Los datos son idénticos (STAR counts del GDC).
# -----------------------------------------------------------------------------

#' Localiza gdc-client (o lo instala en gdc_dir). Devuelve la ruta al ejecutable.
.ensure_gdc_client <- function(gdc_dir) {
  pat <- if (.Platform$OS.type == "windows") "^gdc-client\\.exe$" else "^gdc-client$"
  exe <- list.files(gdc_dir, pattern = pat, recursive = TRUE, full.names = TRUE)
  if (length(exe)) return(exe[1])
  sys <- Sys.info()[["sysname"]]
  url <- switch(sys,
    Windows = "https://gdc.cancer.gov/system/files/public/file/gdc-client_2.3_Windows_x64.zip",
    Darwin  = "https://gdc.cancer.gov/system/files/public/file/gdc-client_2.3_OSX_x64.zip",
    "https://gdc.cancer.gov/system/files/public/file/gdc-client_2.3_Ubuntu_x64.zip")
  zip <- file.path(gdc_dir, "gdc-client.zip")
  utils::download.file(url, zip, mode = "wb", quiet = TRUE)
  utils::unzip(zip, exdir = gdc_dir)
  # algún release anida otro zip
  for (z in list.files(gdc_dir, pattern = "gdc-client.*\\.zip$", full.names = TRUE)) {
    try(utils::unzip(z, exdir = gdc_dir), silent = TRUE)
  }
  exe <- list.files(gdc_dir, pattern = pat, recursive = TRUE, full.names = TRUE)
  if (.Platform$OS.type != "windows" && length(exe)) Sys.chmod(exe[1], "0755")
  stopifnot(length(exe) >= 1)
  exe[1]
}

#' Lee un .tsv STAR augmented del GDC: descarta filas de estadística (N_*) y
#' devuelve gene_id/gene_name/gene_type/unstranded (conteos crudos).
.read_star_counts <- function(f) {
  dt <- data.table::fread(f, sep = "\t", header = TRUE, data.table = FALSE,
                          showProgress = FALSE)
  dt[!grepl("^N_", dt$gene_id),
     c("gene_id", "gene_name", "gene_type", "unstranded"), drop = FALSE]
}

#' Descarga TCGA-LUAD STAR counts con gdc-client y construye el SE directamente.
.gdc_download_prepare <- function(params) {
  old_to <- getOption("timeout"); on.exit(options(timeout = old_to), add = TRUE)
  options(timeout = 3600)
  local_base <- Sys.getenv("LOCALAPPDATA")
  if (local_base == "") local_base <- tempdir()
  gdc_dir  <- params$run$gdc_download_dir %||% file.path(local_base, "luad_gdc")
  data_dir <- file.path(gdc_dir, "files")
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

  # --- catálogo de archivos (metadata; GDCquery sí funciona) ------------------
  query <- TCGAbiolinks::GDCquery(
    project       = params$acquisition$tcga_project,
    data.category = "Transcriptome Profiling",
    data.type     = "Gene Expression Quantification",
    workflow.type = params$acquisition$gdc_workflow)
  res <- TCGAbiolinks::getResults(query)
  # tumor primario + 1 archivo por paciente (determinista)
  res <- res[substr(res$cases, 14, 15) %in% params$acquisition$sample_type_codes, ]
  res <- res[order(res$cases), ]
  res <- res[!duplicated(substr(res$cases, 1, 12)), ]
  if (isTRUE(params$run$smoke_test)) {
    res <- res[seq_len(min(nrow(res), 2L * as.integer(params$run$smoke_n))), ]
  }
  msg(sprintf("TCGA: %d archivos STAR a descargar (gdc-client)", nrow(res)))

  # --- manifiesto y descarga con gdc-client -----------------------------------
  mani <- data.frame(id = res$id, filename = res$file_name, md5 = res$md5sum,
                     size = res$file_size, state = res$state,
                     stringsAsFactors = FALSE)
  mf <- file.path(gdc_dir, "manifest.txt")
  utils::write.table(mani, mf, sep = "\t", quote = FALSE, row.names = FALSE)
  exe <- .ensure_gdc_client(gdc_dir)
  system2(exe, c("download", "-m", mf, "-d", data_dir,
                 "--retry-amount", "5", "-n", "6"),
          stdout = FALSE, stderr = FALSE)

  # --- construir el SummarizedExperiment leyendo los .tsv ---------------------
  files <- file.path(data_dir, res$id, res$file_name)
  ok <- file.exists(files)
  if (any(!ok)) msg(sprintf("aviso: faltan %d/%d archivos tras la descarga",
                            sum(!ok), length(ok)))
  res <- res[ok, ]; files <- files[ok]
  stopifnot(length(files) > 0)

  first  <- .read_star_counts(files[1])
  genes  <- first$gene_id
  counts <- matrix(0L, nrow = length(genes), ncol = length(files),
                   dimnames = list(genes, res$cases))
  for (i in seq_along(files)) {
    dt <- .read_star_counts(files[i])
    counts[, i] <- as.integer(dt$unstranded[match(genes, dt$gene_id)])
  }
  rd <- S4Vectors::DataFrame(gene_id = first$gene_id, gene_name = first$gene_name,
                             gene_type = first$gene_type)
  se <- SummarizedExperiment::SummarizedExperiment(
    assays  = list(unstranded = counts),
    rowData = rd,
    colData = S4Vectors::DataFrame(barcode = res$cases))
  rownames(se) <- genes
  se
}

download_tcga_luad <- function(params) {
  # cache con sufijo por modo (smoke ~2*smoke_n muestras; full ~585).
  suffix <- if (isTRUE(params$run$smoke_test)) "_smoke" else ""
  cache <- normalizePath(
    file.path(params$run$cache_dir, paste0("tcga_luad_se", suffix, ".rds")),
    mustWork = FALSE)
  cache_rds(.gdc_download_prepare(params), cache)   # expr perezoso: no corre si HIT
}

#' Extrae la matriz de conteos (assay 'unstranded' de STAR) genes x muestras.
#' Devuelve conteos crudos (el filtrado/normalización va en 03_qc.R).
tcga_counts <- function(se) {
  assay_name <- if ("unstranded" %in% SummarizedExperiment::assayNames(se))
    "unstranded" else SummarizedExperiment::assayNames(se)[1]
  m <- SummarizedExperiment::assay(se, assay_name)
  storage.mode(m) <- "integer"
  m
}

#' Mapa Ensembl -> símbolo desde el rowData del SE (columna gene_name),
#' devuelto como data.frame (ensembl, symbol, gene_type). El ID de fila puede
#' venir con sufijo de versión (ENSG...\.NN); se conserva tal cual y también
#' su versión sin sufijo para empatar con firmas basadas en Ensembl.
tcga_gene_map <- function(se) {
  rd <- as.data.frame(SummarizedExperiment::rowData(se))
  data.frame(
    ensembl      = rownames(se),
    ensembl_base = sub("\\..*$", "", rownames(se)),
    symbol       = rd$gene_name %||% rd$external_gene_name,
    gene_type    = rd$gene_type,
    stringsAsFactors = FALSE
  )
}

# --- TCGA-CDR (endpoint clínico curado) -------------------------------------

#' Lee TCGA-CDR (Liu et al., Cell 2018). Intenta descargar el .xlsx si no existe.
#' Si la descarga falla, detiene con instrucciones de descarga manual (el UUID
#' del GDC puede cambiar; el objetivo es no continuar en silencio con datos malos).
read_tcga_cdr <- function(params) {
  cfg  <- params$acquisition$cdr
  dest <- file.path(params$run$cache_dir, cfg$file)
  if (!file.exists(dest)) {
    dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
    ok <- tryCatch({
      utils::download.file(cfg$url, dest, mode = "wb", quiet = TRUE); TRUE
    }, error = function(e) FALSE)
    if (!isTRUE(ok) || !file.exists(dest) || file.info(dest)$size < 1e4) {
      stop(sprintf(paste0(
        "No se pudo descargar TCGA-CDR automaticamente.\n",
        "Descarga 'TCGA-CDR-SupplementalTableS1.xlsx' desde la pagina PanCanAtlas:\n",
        "  https://gdc.cancer.gov/about-data/publications/pancanatlas\n",
        "y coloca el archivo en: %s"), dest))
    }
  }
  na_strings <- c("", "NA", "[Not Available]", "[Not Applicable]",
                  "[Discrepancy]", "[Unknown]", "#N/A")
  cdr <- readxl::read_excel(dest, sheet = cfg$sheet, na = na_strings)
  cdr <- as.data.frame(cdr)
  cdr[cdr$type == "LUAD", , drop = FALSE]
}

# --- Ensamblado de la tabla clínica TCGA ------------------------------------

#' Construye la tabla clínica tidy a nivel de muestra: 1 tumor primario por
#' paciente, unido al endpoint curado. Deriva estadio agrupado y variables de
#' lote (TSS, placa) del barcode. Filtra tiempos no positivos o estado ausente
#' (fiabilidad del endpoint). NO estandariza edad aqui (se hace solo en training).
build_clinical <- function(se, cdr, params) {
  endpoint <- params$acquisition$endpoint      # "OS" | "PFI"
  time_col <- paste0(endpoint, ".time")
  ev_col   <- endpoint
  stopifnot(all(c(time_col, ev_col, "bcr_patient_barcode") %in% names(cdr)))

  bc <- colnames(se)
  meta <- data.frame(
    sample_id   = bc,
    patient     = substr(bc, 1, 12),
    sample_type = substr(bc, 14, 15),
    tss         = substr(bc, 6, 7),
    plate       = substr(bc, 22, 25),
    stringsAsFactors = FALSE
  )
  meta <- meta[meta$sample_type %in% params$acquisition$sample_type_codes, ]
  meta <- meta[order(meta$sample_id), ]
  meta <- meta[!duplicated(meta$patient), ]     # 1 muestra/paciente (determinista)

  clin <- data.frame(
    patient   = cdr$bcr_patient_barcode,
    age       = suppressWarnings(as.numeric(cdr$age_at_initial_pathologic_diagnosis)),
    gender    = tolower(as.character(cdr$gender)),
    stage_raw = as.character(cdr$ajcc_pathologic_tumor_stage),
    time      = suppressWarnings(as.numeric(cdr[[time_col]])),
    status    = suppressWarnings(as.integer(cdr[[ev_col]])),
    stringsAsFactors = FALSE
  )

  df <- merge(meta, clin, by = "patient")
  df$stage_group <- collapse_stage(df$stage_raw)
  df$gender <- factor(df$gender, levels = c("female", "male"))
  df$endpoint <- endpoint

  keep <- is.finite(df$time) & df$time > 0 & df$status %in% c(0L, 1L)
  msg(sprintf("clinica TCGA: %d muestras con endpoint %s valido (de %d)",
              sum(keep), endpoint, nrow(df)))
  df[keep, , drop = FALSE]
}

# --- Cohorte externa GEO -----------------------------------------------------

#' Descarga el ExpressionSet de GEO (cacheado). getGPL=TRUE trae la anotacion
#' de plataforma para mapear sondas -> simbolo.
download_geo <- function(params) {
  gse <- params$acquisition$geo$gse_id
  cache <- file.path(params$run$cache_dir, paste0(gse, "_eset.rds"))
  cache_rds({
    gl <- GEOquery::getGEO(gse, GSEMatrix = TRUE, getGPL = TRUE,
                           destdir = params$run$cache_dir)
    # Si hay varias plataformas, escoger la declarada en params.
    idx <- which(vapply(gl, function(e) Biobase::annotation(e),
                        character(1)) == params$acquisition$geo$platform)
    if (length(idx) == 0) idx <- 1
    gl[[idx[1]]]
  }, cache)
}

#' Busca en el fenotipo de GEO (pData) la primera columna que empate un patron.
#' Devuelve NULL si no encuentra (el llamador decide como fallar).
.pheno_find <- function(pd, patterns) {
  cn <- names(pd)
  for (p in patterns) {
    hit <- grep(p, cn, ignore.case = TRUE, value = TRUE)
    if (length(hit)) return(pd[[hit[1]]])
  }
  NULL
}

#' Construye la cohorte GEO: expresion (simbolos) + supervivencia armonizada.
#' El parseo del fenotipo depende del GSE; aqui se cubre GSE68465/GSE31210 con
#' patrones flexibles. Si no localiza tiempo o estado, DETIENE pidiendo revisar
#' colnames(Biobase::pData(eset)) (evita inventar un endpoint).
build_geo_cohort <- function(eset, params) {
  pd <- Biobase::pData(eset)

  # --- supervivencia (tiempo en dias; convertir si viene en meses) ----------
  # OJO: para OS hay que usar el tiempo a MUERTE/ultimo contacto, NO el tiempo a
  # progresion. Se priorizan patrones de OS y se EXCLUYE progresion/relapse.
  t_raw <- .pheno_find(pd, c("last_contact_or_death", "months_to_last",
                             "last.*death", "overall.*surv.*month",
                             "os.*month", "survival.*month"))
  t_days <- .pheno_find(pd, c("survival.*day", "os.*day", "time.*day"))
  vital  <- .pheno_find(pd, c("vital.?status", "death", "os_event",
                              "event", "status"))
  if (is.null(vital) || (is.null(t_raw) && is.null(t_days))) {
    stop(paste0(
      "No pude localizar tiempo/estado de supervivencia en el fenotipo GEO.\n",
      "Revisa: colnames(Biobase::pData(eset)) y ajusta build_geo_cohort()."))
  }
  time <- if (!is.null(t_days)) suppressWarnings(as.numeric(as.character(t_days)))
          else suppressWarnings(as.numeric(as.character(t_raw))) * 30.44  # meses->dias

  v <- tolower(as.character(vital))
  status <- ifelse(grepl("dead|decease|1|yes|event", v), 1L,
            ifelse(grepl("alive|living|0|no|censor", v), 0L, NA_integer_))

  age    <- suppressWarnings(as.numeric(as.character(
              .pheno_find(pd, c("age")))))
  gender <- tolower(as.character(.pheno_find(pd, c("sex", "gender"))))
  gender <- ifelse(grepl("^f", gender), "female",
             ifelse(grepl("^m", gender), "male", NA))
  stage_raw <- as.character(.pheno_find(pd, c("stage", "tnm")))
  # Intentar AJCC directo (p. ej. "Stage IB"); si no aplica y parece TNM
  # ("pN0pT1"), derivar el grupo desde TNM (aprox. AJCC 7ª, M0).
  stage_group <- collapse_stage(stage_raw)
  if (all(is.na(stage_group)) && !is.null(stage_raw) &&
      any(grepl("p[TN]", stage_raw, ignore.case = TRUE))) {
    stage_group <- tnm_to_stage_group(stage_raw)
    msg("GEO: estadio derivado de TNM (aprox. AJCC 7ª, M0)")
  }

  pheno <- data.frame(
    sample_id   = rownames(pd),
    time        = time,
    status      = status,
    age         = age,
    gender      = factor(gender, levels = c("female", "male")),
    stage_group = stage_group,
    stringsAsFactors = FALSE
  )
  keep <- is.finite(pheno$time) & pheno$time > 0 & pheno$status %in% c(0L, 1L)
  pheno <- pheno[keep, , drop = FALSE]

  # --- expresion mapeada a simbolo (colapsando sondas por maxima varianza) ---
  expr <- Biobase::exprs(eset)[, pheno$sample_id, drop = FALSE]
  fd <- Biobase::fData(eset)
  sym_col <- intersect(c("Gene Symbol", "Gene symbol", "GENE_SYMBOL",
                         "Symbol", "gene_symbol"), names(fd))
  if (length(sym_col) == 0) {
    stop(paste0("No hay columna de simbolo en fData(eset). ",
                "Revisa colnames(Biobase::fData(eset))."))
  }
  symbols <- as.character(fd[[sym_col[1]]])
  # Los campos multi-gen tipo 'A /// B' se resuelven al primer simbolo.
  symbols <- sub("\\s*///.*$", "", symbols)
  expr_sym <- collapse_probes_to_symbols(expr, symbols)

  list(expr = expr_sym, pheno = pheno, gse = params$acquisition$geo$gse_id)
}

#' Colapsa una matriz sonda x muestra a simbolo x muestra quedandose, por
#' simbolo, con la sonda de mayor varianza (criterio estandar tipo collapseRows).
collapse_probes_to_symbols <- function(expr, symbols) {
  keep <- !is.na(symbols) & symbols != "" & symbols != "---"
  expr <- expr[keep, , drop = FALSE]
  symbols <- symbols[keep]
  v <- matrixStats::rowVars(as.matrix(expr))
  ord <- order(v, decreasing = TRUE)
  expr <- expr[ord, , drop = FALSE]
  symbols <- symbols[ord]
  first <- !duplicated(symbols)
  out <- expr[first, , drop = FALSE]
  rownames(out) <- symbols[first]
  out
}
