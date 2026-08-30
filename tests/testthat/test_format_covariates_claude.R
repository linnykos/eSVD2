# UNIT_TEST_PLAN.md section 2.1 -- T-FMT-01 .. T-FMT-09.
#
# `format_covariates()` has zero coverage today, and every downstream matrix
# depends on its column names and ordering. Two of the tests below are expected
# to FAIL against the current code: they pin documentation/code mismatches that
# section 2.1 records as defects, and no code has been changed yet.

.fmt_dat <- function(n = 12, p = 5){
  set.seed(10)
  matrix(stats::rpois(n * p, lambda = 5) + 1L, nrow = n, ncol = p,
         dimnames = list(paste0("cell_", seq_len(n)), paste0("gene_", seq_len(p))))
}

test_that("T-FMT-01: column order is Intercept, Log_UMI, numerics, then factors", {
  dat <- .fmt_dat()
  covariate_df <- data.frame(Age = stats::rnorm(nrow(dat)),
                             Group = factor(rep(c("a", "b", "c"), length.out = nrow(dat))))

  res <- format_covariates(dat = dat, covariate_df = covariate_df)

  expect_equal(colnames(res)[1], "Intercept")
  expect_equal(colnames(res)[2], "Log_UMI")
  # The factor indicators come last, in `levels()` order with the first dropped.
  expect_equal(utils::tail(colnames(res), 2), c("Group_b", "Group_c"))
  expect_equal(nrow(res), nrow(dat))
})

## Section 2.1 defect: the roxygen says "dropping the last level"; the code
## drops the FIRST. One of the two must change. This test asserts the behaviour
## that agrees with `stats::model.matrix`, which is the convention every reader
## will assume, so it is expected to PASS and the documentation is what moves.
test_that("T-FMT-02: a 3-level factor yields 2 indicators, first level dropped", {
  dat <- .fmt_dat(n = 6)
  group_vec <- factor(c("a", "a", "b", "b", "c", "c"))
  covariate_df <- data.frame(Group = group_vec)

  res <- format_covariates(dat = dat, covariate_df = covariate_df)
  indicator_mat <- res[, c("Group_b", "Group_c"), drop = FALSE]

  model_mat <- stats::model.matrix(~ group_vec)[, -1, drop = FALSE]
  expect_equal(unname(indicator_mat), unname(model_mat))
  expect_equal(sum(grepl("^Group_", colnames(res))), 2)
})

test_that("T-FMT-03: variables_enumerate_all keeps every level", {
  dat <- .fmt_dat(n = 6)
  covariate_df <- data.frame(Group = factor(c("a", "a", "b", "b", "c", "c")))

  res <- format_covariates(dat = dat,
                           covariate_df = covariate_df,
                           variables_enumerate_all = "Group")

  expect_equal(sum(grepl("^Group_", colnames(res))), 3)
  expect_true(all(c("Group_a", "Group_b", "Group_c") %in% colnames(res)))
  # Keeping all levels makes the indicator block sum to 1 in every row, which
  # is collinear with the intercept. That is the caller's problem, but the
  # test records that the branch does what it says.
  expect_equal(unname(rowSums(res[, c("Group_a", "Group_b", "Group_c")])),
               rep(1, nrow(dat)))
})

test_that("T-FMT-04: Log_UMI equals log(rowSums(dat)) for dense and sparse", {
  dat <- .fmt_dat()
  covariate_df <- data.frame(Age = stats::rnorm(nrow(dat)))

  res_dense <- format_covariates(dat = dat, covariate_df = covariate_df)
  expect_equal(unname(res_dense[, "Log_UMI"]),
               unname(log(Matrix::rowSums(dat))))

  dat_sparse <- methods::as(methods::as(dat * 1.0, "dMatrix"), "CsparseMatrix")
  res_sparse <- format_covariates(dat = dat_sparse, covariate_df = covariate_df)
  expect_equal(unname(res_sparse[, "Log_UMI"]),
               unname(log(Matrix::rowSums(dat_sparse))))
  expect_equal(res_dense, res_sparse)
})

## NEW FINDING, not in UNIT_TEST_PLAN.md section 2.1 and not in
## CRAN_READINESS.md: "rescale" does NOT standardize to unit variance.
##
## `scale(x, center = FALSE, scale = TRUE)` divides by the root-mean-square
## sqrt(sum(x^2)/(n-1)), not by the standard deviation. Those coincide only when
## the column is already centered. For `Age` with mean 40 and sd 8 the RMS is
## about 40.8, so the rescaled column has sd about 0.2 -- not 1.
##
## The consequence is that the scaled magnitude of a covariate depends on its
## MEAN, so `Age` in years and `Age` in months rescale identically (both are
## divided by their own RMS) but `Age` and `Age - 40` do not. Whether that is
## intended is a modelling question; see the report. This test pins the actual
## semantics so a future edit to the `scale()` arguments is caught either way.
test_that("T-FMT-05: a rescaled numeric is divided by its RMS, not its sd", {
  dat <- .fmt_dat()
  age_vec <- stats::rnorm(nrow(dat), mean = 40, sd = 8)
  covariate_df <- data.frame(Age = age_vec)

  res <- format_covariates(dat = dat,
                           covariate_df = covariate_df,
                           bool_center = FALSE,
                           rescale_numeric_variables = "Age")

  root_mean_square <- sqrt(sum(age_vec^2) / (length(age_vec) - 1))
  expect_equal(unname(res[, "Age"]), unname(age_vec / root_mean_square),
               tolerance = 1e-10)

  # The paper is explicit about scaling WITHOUT centering: centering here would
  # change the interpretation of the intercept. That much the code gets right.
  expect_false(isTRUE(all.equal(mean(res[, "Age"]), 0)))
  expect_true(all(res[, "Age"] > 0))

  # ... but it is emphatically not unit variance, which is what a reader of the
  # word "rescales" would assume.
  expect_false(isTRUE(all.equal(stats::sd(res[, "Age"]), 1,
                                tolerance = 1e-3)))
})

test_that("T-FMT-05b: rescaling is equivariant to the units of the covariate", {
  dat <- .fmt_dat()
  age_vec <- stats::rnorm(nrow(dat), mean = 40, sd = 8)

  res_years <- format_covariates(dat = dat,
                                 covariate_df = data.frame(Age = age_vec),
                                 rescale_numeric_variables = "Age")
  res_months <- format_covariates(dat = dat,
                                  covariate_df = data.frame(Age = 12 * age_vec),
                                  rescale_numeric_variables = "Age")

  # This is the property the RMS division actually delivers, and it is a real
  # one: the fitted coefficient does not depend on the unit the user recorded.
  expect_equal(res_years[, "Age"], res_months[, "Age"], tolerance = 1e-10)
})

## Section 2.1 defect: the roxygen says the function "rescales all the numerical
## variables"; the code rescales only those named in
## `rescale_numeric_variables`. This test asserts the CODE's behaviour, since it
## is the more useful of the two, so it should pass and the docs should move.
test_that("T-FMT-06: a numeric not named in rescale_numeric_variables is untouched", {
  dat <- .fmt_dat()
  age_vec <- stats::rnorm(nrow(dat), mean = 40, sd = 8)
  weight_vec <- stats::rnorm(nrow(dat), mean = 70, sd = 12)
  covariate_df <- data.frame(Age = age_vec, Weight = weight_vec)

  res <- format_covariates(dat = dat,
                           covariate_df = covariate_df,
                           rescale_numeric_variables = "Age")

  expect_equal(unname(res[, "Weight"]), unname(weight_vec))
  expect_false(isTRUE(all.equal(stats::sd(res[, "Weight"]), 1)))
})

test_that("T-FMT-07: a row-count mismatch errors informatively", {
  dat <- .fmt_dat(n = 12)
  covariate_df <- data.frame(Age = stats::rnorm(11))

  expect_error(format_covariates(dat = dat, covariate_df = covariate_df),
               regexp = "nrow")
})

test_that("T-FMT-08: a single-level factor errors, naming the variable", {
  dat <- .fmt_dat(n = 6)
  covariate_df <- data.frame(Group = factor(rep("a", 6)))

  # `eSVD()` pre-filters these, so the error is user-facing only for direct
  # callers -- but it should still say which variable is at fault.
  expect_error(format_covariates(dat = dat, covariate_df = covariate_df),
               regexp = "Group")
})

test_that("T-FMT-09: a factor level with spaces and parentheses survives verbatim", {
  dat <- .fmt_dat(n = 6)
  covariate_df <- data.frame(
    Diagnosis = factor(c("control", "control", "ASD (severe)", "ASD (severe)",
                         "control", "ASD (severe)"))
  )

  res <- format_covariates(dat = dat, covariate_df = covariate_df)

  # This is the input that breaks `reparameterization_esvd_covariates()`
  # downstream (CRAN_READINESS.md section 1.3), because `as.data.frame()` there
  # applies `make.names()`. Asserting the name here localises the bug to the
  # right function.
  expect_true("Diagnosis_control" %in% colnames(res))
  expect_false(any(grepl("Diagnosis.control", colnames(res), fixed = TRUE)))
})
