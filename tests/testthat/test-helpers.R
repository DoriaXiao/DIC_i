test_that("print.DICi prints all components and returns input invisibly", {
  set.seed(8)
  ll <- matrix(rnorm(100 * 5, -2, 0.5), 100, 5)
  res <- compute_dic_i(ll)

  expect_output(print(res), "DIC_i")
  expect_output(print(res), "p_V")
  expect_output(print(res), "E\\[D\\]")
  expect_output(print(res, digits = 4),
                sprintf("DIC_i  = %.4f", res$dic_i), fixed = TRUE)
  expect_invisible(print(res))
  expect_identical(withVisible(print(res))$value, res)
})

test_that("dici_example lists the bundled files", {
  expect_output(files <- dici_example(), "Available example files")
  expect_setequal(files, c("example_fa.R", "fa_marginal.stan"))
})

test_that("dici_example(copy = FALSE) returns an existing path", {
  path <- dici_example("fa_marginal.stan", copy = FALSE)
  expect_true(file.exists(path))
  expect_equal(basename(path), "fa_marginal.stan")
})

test_that("dici_example(copy = TRUE) copies into the working directory", {
  tmp <- tempfile("dici_example_")
  dir.create(tmp)
  old <- setwd(tmp)
  on.exit({
    setwd(old)
    unlink(tmp, recursive = TRUE)
  }, add = TRUE)

  expect_message(dst <- dici_example("example_fa.R"), "Copied")
  expect_true(file.exists(file.path(tmp, "example_fa.R")))
  expect_equal(normalizePath(dst), normalizePath(file.path(tmp, "example_fa.R")))

  # Second copy overwrites and says so
  expect_message(dici_example("example_fa.R"), "Overwriting")
})

test_that("dici_example errors on unknown file", {
  expect_error(dici_example("no_such_file.stan", copy = FALSE), "not found")
})

test_that("dic_i_from_cmdstanr errors informatively without cmdstanr", {
  skip_if(requireNamespace("cmdstanr", quietly = TRUE),
          "cmdstanr is installed")
  expect_error(dic_i_from_cmdstanr(list()), "cmdstanr")
})

test_that("dic_i_from_cmdstanr matches compute_dic_i on a mock fit", {
  skip_if_not_installed("cmdstanr")
  skip_if_not_installed("posterior")

  set.seed(9)
  ll <- matrix(rnorm(200 * 6, -2, 0.5), 200, 6,
               dimnames = list(NULL, sprintf("log_lik[%d]", 1:6)))
  # Minimal stand-in for a CmdStanMCMC object: only $draws() is used
  mock_fit <- list(
    draws = function(variables, format) posterior::as_draws_matrix(ll)
  )

  res <- dic_i_from_cmdstanr(mock_fit)
  ref <- compute_dic_i(unname(ll))

  expect_s3_class(res, "DICi")
  expect_equal(res$dic_i, ref$dic_i)
  expect_equal(res$p_v,   ref$p_v)
  expect_equal(res$e_d,   ref$e_d)
})
