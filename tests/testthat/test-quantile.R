## Quantile tables are no longer stored in the fits: the fields mod$pred_q, ...
## rebuild them on access from the stored link-scale mean/SD and the holdout
## calibration, at the 15 default levels; predict() takes 'probs'.

test_that("pred_q is rebuilt on access, and predict() gives other levels", {
  d  <- sim_spatial(n = 120); i0 <- 1:15
  hv <- quiet(cf_lm_hv(y = d$y, x = d$x, coords = d$coords))
  m  <- quiet(cf_lm(y = d$y, x = d$x, coords = d$coords, x0 = d$x[i0, ],
                    coords0 = d$coords[i0, ] + 0.01, mod_hv = hv))
  expect_null(.subset2(m, "pred_q"))                       # not stored
  q  <- m$pred_q
  expect_identical(m[["pred_q"]], q)
  expect_equal(dim(q), c(120L, 15L))
  expect_equal(nrow(m$pred0_q), 15L)
  expect_false("probs" %in% names(formals(cf_lm)))
  ## other levels through predict()
  q2 <- predict(m, probs = c(0.025, 0.5, 0.975))[, -(1:2)]
  expect_named(q2, c("q0.025", "q0.5", "q0.975"))
  expect_identical(q2$q0.025, q$q0.025)
  ## the 95% interval quantiles are consistent with pred / pred_sd (Gaussian)
  expect_equal(q2$q0.5, m$pred$pred)
  expect_equal(q2$q0.975 - q2$q0.5, qnorm(0.975) * m$pred$pred_sd)
  ## the signal quantiles describe the mean, not a new observation
  qs <- predict(m, probs = c(0.025, 0.975), se_type = "mean")
  expect_equal((qs$q0.975 - qs$q0.025) / 2, qnorm(0.975) * m$pred_signal$pred_sd)
  ## other fields keep the default `$` behaviour, including partial matching
  expect_identical(m$beta, .subset2(m, "beta"))
  expect_identical(m$bet, .subset2(m, "beta"))
  expect_identical(m$pred_q_signal[, c("q0.025", "q0.975")], qs[, c("q0.025", "q0.975")])
  expect_error(predict(m, probs = 1.2), "between 0 and 1")
})

test_that("pred_q_signal is NULL and quantiles are Gaussian for se_type = 'mean'", {
  d  <- sim_spatial(n = 120)
  hv <- quiet(cf_lm_hv(y = d$y, x = d$x, coords = d$coords))
  m  <- quiet(cf_lm(y = d$y, x = d$x, coords = d$coords, mod_hv = hv, se_type = "mean"))
  expect_null(m$pred_q_signal)
  expect_null(m$pred0_q)                                   # no prediction sites
})

test_that("count families rebuild the calibrated observation quantiles", {
  d  <- sim_spatial(n = 120); set.seed(4)
  y  <- rpois(120, exp(0.3 * d$x$v1 + 0.6 * d$field))
  hv <- quiet(cf_glm_hv(y = y, x = d$x, coords = d$coords, family = poisson()))
  m  <- quiet(cf_glm(y = y, x = d$x, coords = d$coords, mod_hv = hv))
  q  <- predict(m, probs = c(0.025, 0.5, 0.975))[, c("q0.025", "q0.5", "q0.975")]
  expect_true(all(q == round(q)))                          # integer-valued (negative binomial)
  expect_true(all(q$q0.025 <= q$q0.5 & q$q0.5 <= q$q0.975))
})

test_that("fits saved by spCF <= 0.2.1 return their stored quantiles", {
  d  <- sim_spatial(n = 120)
  hv <- quiet(cf_lm_hv(y = d$y, x = d$x, coords = d$coords))
  m  <- quiet(cf_lm(y = d$y, x = d$x, coords = d$coords, mod_hv = hv))
  q  <- m$pred_q
  old <- unclass(m); old$other$qspec <- NULL; old$other$pcore <- NULL; old$pred_q <- q; class(old) <- "cf_lm"
  expect_identical(old$pred_q, q)
  expect_identical(predict(old, probs = c(0.025, 0.975))[, c("q0.025", "q0.975")], q[, c("q0.025", "q0.975")])
  expect_error(predict(old, probs = 0.33), "only at the levels")
  expect_error(predict(old, coords0 = d$coords[1:3, ]), "refit")
})

test_that("keep_scales = FALSE drops the scale-wise processes only", {
  d  <- sim_spatial(n = 120); i0 <- 1:10
  hv <- quiet(cf_lm_hv(y = d$y, x = d$x, coords = d$coords))
  a  <- quiet(cf_lm(y = d$y, x = d$x, coords = d$coords, x0 = d$x[i0, ], coords0 = d$coords[i0, ], mod_hv = hv))
  b  <- quiet(cf_lm(y = d$y, x = d$x, coords = d$coords, x0 = d$x[i0, ], coords0 = d$coords[i0, ], mod_hv = hv,
                    keep_scales = FALSE))
  expect_false(is.null(a$Z))
  expect_null(b$Z); expect_null(b$Z_sd); expect_null(b$Z0); expect_null(b$Z0_sd)
  expect_identical(a$pred, b$pred); expect_identical(a$pred0, b$pred0)
  expect_identical(a$beta, b$beta); expect_identical(a$sd_summary, b$sd_summary)
  expect_identical(a$pred_q, b$pred_q)
  expect_error(sp_scalewise(b), "keep_scales")
  expect_lt(as.numeric(object.size(b)), as.numeric(object.size(a)))

  s  <- sim_spacetime(ns = 40, nt = 4)
  hv <- quiet(cf_dglm_hv(y = s$y, x = s$x, coords = s$coords, time = s$time))
  a  <- quiet(cf_dglm(y = s$y, x = s$x, coords = s$coords, time = s$time, mod_hv = hv))
  b  <- quiet(cf_dglm(y = s$y, x = s$x, coords = s$coords, time = s$time, mod_hv = hv, keep_scales = FALSE))
  expect_null(b$Z); expect_identical(a$pred, b$pred); expect_identical(a$pred_q, b$pred_q)
})

test_that("prediction tables carry automatic (integer) row names", {
  ## named prediction vectors used to become character row names, which made
  ## a two-column table about five times larger than its numbers
  d  <- sim_spatial(n = 120); set.seed(4)
  y  <- rpois(120, exp(0.3 * d$x$v1 + 0.6 * d$field))
  hv <- quiet(cf_glm_hv(y = y, x = d$x, coords = d$coords, family = poisson()))
  m  <- quiet(cf_glm(y = y, x = d$x, coords = d$coords, x0 = d$x[1:5, ], coords0 = d$coords[1:5, ], mod_hv = hv))
  expect_lt(.row_names_info(m$pred), 0)
  expect_lt(.row_names_info(m$pred0), 0)
  s  <- sim_spacetime(ns = 40, nt = 4)
  hv <- quiet(cf_dglm_hv(y = s$y, x = s$x, coords = s$coords, time = s$time))
  m  <- quiet(cf_dglm(y = s$y, x = s$x, coords = s$coords, time = s$time, mod_hv = hv))
  expect_lt(.row_names_info(m$pred), 0)
})
