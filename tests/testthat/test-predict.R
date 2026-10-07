## predict() at new sites uses only the knot states stored in the fit, and must
## agree exactly with fitting the model with the same sites as coords0.

## Predictions computed in different batches (or on a refit, or after saveRDS)
## go through matrix products of different sizes, which alternative BLAS
## libraries (BLIS, OpenBLAS) evaluate in a different order; the results then
## agree only up to rounding (~1e-16), so compare with a small tolerance.
expect_same <- function(object, expected, ...)
  expect_equal(object, expected, tolerance = 1e-10, ...)

test_that("predict.cf_lm agrees with cf_lm(coords0 = ...) and is batch-invariant", {
  d  <- sim_spatial(n = 150); set.seed(8)
  g  <- data.frame(px = runif(40), py = runif(40)); x0 <- data.frame(v1 = rnorm(40), v2 = runif(40))
  hv <- quiet(cf_lm_hv(y = d$y, x = d$x, coords = d$coords))
  mw <- quiet(cf_lm(y = d$y, x = d$x, coords = d$coords, x0 = x0, coords0 = g, mod_hv = hv))
  mo <- quiet(cf_lm(y = d$y, x = d$x, coords = d$coords, mod_hv = hv))
  p  <- predict(mo, x0 = x0, coords0 = g)
  expect_same(as.matrix(p), as.matrix(cbind(mw$pred0, mw$pred0_q)))
  pa <- rbind(predict(mo, x0 = x0[1:15, ], coords0 = g[1:15, ]),
              predict(mo, x0 = x0[16:40, ], coords0 = g[16:40, ]))
  expect_same(as.matrix(pa), as.matrix(p))
  ## levels, the mean, and the sample sites
  expect_named(predict(mo, x0 = x0, coords0 = g, probs = c(0.1, 0.9)), c("pred", "pred_sd", "q0.1", "q0.9"))
  pm <- predict(mo, x0 = x0, coords0 = g, se_type = "mean")
  expect_same(pm$pred, p$pred)
  expect_same(pm$pred_sd, mw$pred0_signal$pred_sd)
  expect_same(predict(mo)$pred, mo$pred$pred)
  ## a saved and reloaded fit predicts the same
  f <- tempfile(fileext = ".rds"); saveRDS(mo, f)
  expect_same(predict(readRDS(f), x0 = x0, coords0 = g), p)
})

test_that("predict.cf_glm agrees with cf_glm(coords0 = ...), with an offset", {
  d   <- sim_spatial(n = 150); set.seed(9)
  off <- log(runif(150, 1, 3)); y <- rpois(150, exp(off + 0.3 * d$x$v1 + 0.6 * d$field))
  g   <- data.frame(px = runif(30), py = runif(30)); x0 <- data.frame(v1 = rnorm(30), v2 = runif(30))
  off0 <- log(runif(30, 1, 3))
  hv  <- quiet(cf_glm_hv(y = y, x = d$x, coords = d$coords, offset = off, family = poisson()))
  mw  <- quiet(cf_glm(y = y, x = d$x, coords = d$coords, offset = off, x0 = x0, coords0 = g,
                      offset0 = off0, mod_hv = hv))
  mo  <- quiet(cf_glm(y = y, x = d$x, coords = d$coords, offset = off, mod_hv = hv))
  expect_same(as.matrix(predict(mo, x0 = x0, coords0 = g, offset0 = off0)),
                   as.matrix(cbind(mw$pred0, mw$pred0_q)))
})

test_that("predict checks its inputs", {
  d  <- sim_spatial(n = 150)
  hv <- quiet(cf_lm_hv(y = d$y, x = d$x, coords = d$coords))
  mo <- quiet(cf_lm(y = d$y, x = d$x, coords = d$coords, mod_hv = hv, se_type = "mean"))
  g  <- d$coords[1:5, ]
  expect_error(predict(mo, coords0 = g), "'x0' must be provided")
  expect_error(predict(mo, x0 = d$x[1:4, ], coords0 = g), "row")
  expect_error(predict(mo, x0 = d$x[1:5, 1, drop = FALSE], coords0 = g), "column")
  expect_error(predict(mo, x0 = d$x[1:5, ], coords0 = g, se_type = "prediction"), "se_type")
  expect_error(predict(mo, x0 = d$x[1:5, ], coords0 = g, probs = 0), "between 0 and 1")
})

test_that("predict.cf_dglm agrees with cf_dglm(time0 = ...) at training, interior and future times", {
  set.seed(11); s <- sim_spacetime(ns = 40, nt = 6); keep <- s$time != 3     # time 3 not observed
  y <- s$y[keep]; X <- s$x[keep, , drop = FALSE]; C <- s$coords[keep, ]; tt <- s$time[keep]
  hv <- quiet(cf_dglm_hv(y = y, x = X, coords = C, time = tt))
  mo <- quiet(cf_dglm(y = y, x = X, coords = C, time = tt, mod_hv = hv))
  set.seed(2); c0 <- data.frame(px = runif(12), py = runif(12)); x0 <- data.frame(v1 = rnorm(12))
  t0 <- rep(c(0, 2, 3, 3.5, 6, 8), 2)
  mw <- quiet(cf_dglm(y = y, x = X, coords = C, time = tt, mod_hv = hv, x0 = x0, coords0 = c0, time0 = t0))
  expect_same(mo$pred, mw$pred)                                  # the fit does not depend on time0
  p  <- predict(mo, x0 = x0, coords0 = c0, time0 = t0)
  expect_same(as.matrix(p), as.matrix(cbind(mw$pred0, mw$pred0_q)))
  pa <- rbind(predict(mo, x0 = x0[1:5, , drop = FALSE], coords0 = c0[1:5, ], time0 = t0[1:5]),
              predict(mo, x0 = x0[6:12, , drop = FALSE], coords0 = c0[6:12, ], time0 = t0[6:12]))
  expect_same(as.matrix(pa), as.matrix(p))
  ## the bridge is continuous at the training times, and uncertainty grows
  ## between them and into the future (link scale = response scale here)
  site <- s$coords_uni[rep(5, 4), ]; v0 <- data.frame(v1 = rep(0, 4))
  pt <- predict(mo, x0 = v0, coords0 = site, time0 = c(2, 2 + 1e-6, 4 - 1e-6, 4), se_type = "mean")
  expect_lt(abs(pt$pred[1] - pt$pred[2]), 1e-4); expect_lt(abs(pt$pred[3] - pt$pred[4]), 1e-4)
  pm <- predict(mo, x0 = v0[1:3, , drop = FALSE], coords0 = site[1:3, ], time0 = c(2, 3, 4), se_type = "mean")
  expect_gt(pm$pred_sd[2], max(pm$pred_sd[c(1, 3)]))
  pf <- predict(mo, x0 = v0[1:3, , drop = FALSE], coords0 = site[1:3, ], time0 = c(6, 8, 12), se_type = "mean")
  expect_true(all(diff(pf$pred_sd) > 0))
  expect_error(predict(mo, x0 = x0, coords0 = c0), "time0")
})
