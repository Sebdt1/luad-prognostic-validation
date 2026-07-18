# =============================================================================
# test_utils.R — Tests unitarios (testthat) de las utilidades criticas
# -----------------------------------------------------------------------------
# Cubren la logica pura donde un bug corromperia el analisis en silencio:
# colapso de estadio, escalado train-only, determinismo de la particion y
# ausencia de solapamiento train/test. Ejecutar desde la raiz del proyecto:
#   Rscript -e 'testthat::test_file("tests/test_utils.R")'
# =============================================================================

library(testthat)

# Localiza R/ de forma robusta (ejecutable con test_file desde la raiz o desde
# tests/). Se hace ANTES de sourcear, sin depender de utilidades del proyecto.
find_R_dir <- function() {
  for (p in c("R", file.path("..", "R"))) if (dir.exists(p)) return(normalizePath(p))
  stop("No encuentro el directorio R/ (ejecuta desde la raiz del proyecto).")
}
R_dir <- find_R_dir()

# Cargar solo los modulos con funciones puras (sin dependencias de Bioconductor).
source(file.path(R_dir, "00_utils.R"))
source(file.path(R_dir, "01_split.R"))

# Parametros minimos de prueba
params_stub <- list(
  project = list(seed = 1234),
  run     = list(smoke_test = TRUE, smoke_n = 6),
  split   = list(p_train = 0.7, stratify_by = "stage_group")
)

# --- collapse_stage ----------------------------------------------------------
test_that("collapse_stage mapea correctamente y NA lo no informativo", {
  x <- c("Stage IA", "Stage IB", "Stage IIA", "Stage IIIB", "Stage IV",
         "[Not Available]", "Stage X", NA, "")
  out <- collapse_stage(x)
  expect_equal(as.character(out),
               c("I", "I", "II", "III", "IV", NA, NA, NA, NA))
  expect_equal(levels(out), c("I", "II", "III", "IV"))
})

# --- zscore_fit / zscore_apply ----------------------------------------------
test_that("zscore aprende center/scale en un vector y los aplica a otro", {
  train <- c(1, 2, 3, 4, 5)
  fit <- zscore_fit(train)
  expect_equal(fit$center, 3)
  expect_equal(fit$scale, stats::sd(train))
  # aplicar al MISMO vector -> media 0, sd 1
  z <- zscore_apply(train, fit)
  expect_equal(mean(z), 0, tolerance = 1e-8)
  expect_equal(stats::sd(z), 1, tolerance = 1e-8)
  # escala constante -> no divide por 0
  fit0 <- zscore_fit(c(2, 2, 2))
  expect_equal(fit0$scale, 1)
})

test_that("rank_normalize cae en [0,1] y es monotono", {
  x <- c(10, 30, 20, 40)
  r <- rank_normalize(x)
  expect_true(all(r >= 0 & r <= 1))
  expect_equal(order(r), order(x))
})

# --- %||% --------------------------------------------------------------------
test_that("%||% coalesce nulos y vacios", {
  expect_equal(NULL %||% 5, 5)
  expect_equal(character(0) %||% "x", "x")
  expect_equal(3 %||% 5, 3)
})

# --- make_split: determinismo, sin solapamiento, estratificacion ------------
make_fake_clinical <- function(n = 40) {
  set.seed(99)
  data.frame(
    sample_id = sprintf("S%03d", seq_len(n)),
    stage_group = factor(sample(c("I", "II", "III", "IV"), n, replace = TRUE),
                         levels = c("I", "II", "III", "IV")),
    stringsAsFactors = FALSE
  )
}

test_that("make_split es determinista con la misma semilla", {
  cl <- make_fake_clinical()
  s1 <- make_split(cl, params_stub)
  s2 <- make_split(cl, params_stub)
  expect_identical(s1$train, s2$train)
  expect_identical(s1$test, s2$test)
})

test_that("make_split no solapa train/test y cubre toda la cohorte", {
  cl <- make_fake_clinical()
  s <- make_split(cl, params_stub)
  expect_length(intersect(s$train, s$test), 0)
  expect_setequal(c(s$train, s$test), cl$sample_id)
})

test_that("make_split respeta approx la proporcion de training", {
  cl <- make_fake_clinical(200)
  s <- make_split(cl, params_stub)
  expect_gt(length(s$train) / 200, 0.6)
  expect_lt(length(s$train) / 200, 0.8)
})

# --- subset_to ---------------------------------------------------------------
test_that("subset_to filtra data.frames por sample_id y matrices por columna", {
  cl <- make_fake_clinical(10)
  sub <- subset_to(cl, cl$sample_id[1:3])
  expect_equal(nrow(sub), 3)
  m <- matrix(1:20, nrow = 4, dimnames = list(NULL, sprintf("S%03d", 1:5)))
  msub <- subset_to(m, c("S001", "S003"))
  expect_equal(colnames(msub), c("S001", "S003"))
})

# --- maybe_subsample ---------------------------------------------------------
test_that("maybe_subsample respeta el objetivo y es no-op si smoke=FALSE", {
  ids <- sprintf("S%03d", 1:100)
  set.seed(1)
  sub <- maybe_subsample(ids, params_stub)
  expect_lte(length(sub), 6)
  p_off <- modifyList(params_stub, list(run = list(smoke_test = FALSE, smoke_n = 6)))
  expect_length(maybe_subsample(ids, p_off), 100)
})
