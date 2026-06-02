test_that("compute_dic_i returns correct structure", {
  set.seed(1)
  ll <- matrix(rnorm(200 * 10, -2, 0.5), 200, 10)
  res <- compute_dic_i(ll)

  expect_s3_class(res, "DICi")
  expect_named(res, c("dic_i", "p_v", "e_d"))
  expect_true(is.numeric(res$dic_i))
  expect_true(is.numeric(res$p_v))
  expect_true(is.numeric(res$e_d))
})

test_that("p_v is always non-negative", {
  set.seed(2)
  ll <- matrix(rnorm(100 * 5, -3, 1), 100, 5)
  res <- compute_dic_i(ll)
  expect_gte(res$p_v, 0)
})

test_that("log_lik and deviance_draws give identical results", {
  set.seed(3)
  ll <- matrix(rnorm(500 * 20, -2, 0.5), 500, 20)
  dev <- -2 * rowSums(ll)

  res_ll  <- compute_dic_i(log_lik = ll)
  res_dev <- compute_dic_i(deviance_draws = dev)

  expect_equal(res_ll$dic_i, res_dev$dic_i)
  expect_equal(res_ll$p_v,   res_dev$p_v)
  expect_equal(res_ll$e_d,   res_dev$e_d)
})

test_that("DIC_i = E[D] + p_V", {
  set.seed(4)
  ll <- matrix(rnorm(300 * 15, -1.5, 0.3), 300, 15)
  res <- compute_dic_i(ll)
  expect_equal(res$dic_i, res$e_d + res$p_v)
})

test_that("manual computation matches function output", {
  set.seed(5)
  ll <- matrix(rnorm(400 * 8, -2, 0.4), 400, 8)

  d <- -2 * rowSums(ll)
  expected_e_d <- mean(d)
  expected_p_v <- 0.5 * var(d)
  expected_dic <- expected_e_d + expected_p_v

  res <- compute_dic_i(ll)
  expect_equal(res$e_d,   expected_e_d)
  expect_equal(res$p_v,   expected_p_v)
  expect_equal(res$dic_i, expected_dic)
})

test_that("input validation works", {
  expect_error(compute_dic_i(), "Either")
  expect_error(compute_dic_i(log_lik = "not a matrix"), "must be a matrix")
  expect_error(compute_dic_i(log_lik = matrix(1, 1, 5)), "at least 2 rows")
  expect_error(compute_dic_i(deviance_draws = 42), "at least 2 elements")
  expect_error(compute_dic_i(deviance_draws = "bad"), "numeric vector")
})

test_that("compare_dic_i returns sorted table", {
  set.seed(6)
  ll1 <- matrix(rnorm(200 * 10, -2.0, 0.5), 200, 10)
  ll2 <- matrix(rnorm(200 * 10, -2.5, 0.5), 200, 10)

  f1 <- compute_dic_i(ll1)
  f2 <- compute_dic_i(ll2)

  tbl <- compare_dic_i(good = f1, bad = f2)

  expect_s3_class(tbl, "data.frame")
  expect_equal(nrow(tbl), 2)
  expect_equal(tbl$delta_dic_i[1], 0)
  expect_true(tbl$dic_i[1] <= tbl$dic_i[2])
})

test_that("compare_dic_i validates inputs", {
  set.seed(7)
  ll <- matrix(rnorm(100 * 5), 100, 5)
  f1 <- compute_dic_i(ll)

  expect_error(compare_dic_i(a = f1), "At least two")
  expect_error(compare_dic_i(f1, f1), "must be named")
  expect_error(compare_dic_i(a = f1, b = "not_dici"), "not a DICi object")
})
