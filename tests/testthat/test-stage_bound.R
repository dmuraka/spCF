## Per-stage bounded predictive variance (internal_stage_var.R). The bound is always
## on; the internal option spcf.stage_bound = FALSE restores spCF <= 0.2.1 for comparison.

d  <- sim_spatial(n = 150, seed = 11)
hv <- quiet(cf_lm_hv(y = d$y, x = d$x, coords = d$coords))
## prediction sites from inside the data cloud to far outside it
far <- data.frame(px = c(0.5, 0.5, 1.2, 2, 4), py = c(0.5, 0.5, 1.2, 2, 4))
x0  <- data.frame(v1 = 0, v2 = 0.5)[rep(1, nrow(far)), ]

test_that("stage_bound leaves point predictions and coefficients unchanged", {
  a <- quiet(cf_lm(y = d$y, x = d$x, coords = d$coords, x0 = x0, coords0 = far, mod_hv = hv))
  b <- withr::with_options(list(spcf.stage_bound = FALSE), quiet(cf_lm(y = d$y, x = d$x, coords = d$coords, x0 = x0, coords0 = far, mod_hv = hv)))
  expect_equal(a$pred0$pred, b$pred0$pred)
  expect_equal(a$beta_int_summ$coef, b$beta_int_summ$coef)
})

test_that("bounded field variance is finite, grows away from the data and stays below the sill", {
  m <- quiet(cf_lm(y = d$y, x = d$x, coords = d$coords, x0 = x0, coords0 = far, mod_hv = hv))
  Z0sd <- as.matrix(m$Z0_sd); Z <- as.matrix(m$Z)
  fv0 <- rowSums(Z0sd^2)
  expect_true(all(is.finite(fv0)))
  expect_true(all(fv0 <= stats::var(rowSums(Z)) + 1e-8))
  expect_true(all(diff(m$pred0$pred_sd[2:5]) >= -1e-8))   # non-decreasing with distance
  ## far outside the data every stage is at its cap: the field variance equals the sill
  expect_equal(fv0[5], stats::var(rowSums(Z)), tolerance = 1e-6)
})

test_that("stage caps sum to the sill and the tau equation is solved", {
  Z <- cbind(rnorm(50), rnorm(50, sd = .5)); Z <- Z + rnorm(50)   # correlated stages
  sill <- stats::var(rowSums(Z)); caps <- spCF:::.spcf_stage_caps(Z, sill)
  expect_equal(sum(caps), sill)
  P <- matrix(runif(100, .1, 2), 50, 2); num <- 0.5 * sill
  t <- spCF:::.spcf_stage_tau_raw(P, caps, num)
  expect_equal(mean(rowSums(spCF:::.spcf_stage_var(P, caps, t))), num, tolerance = 1e-6)
  expect_equal(spCF:::.spcf_stage_tau_raw(P, caps, 2 * sill), 100)
})

test_that("cf_glm uses the bound for poisson and ignores it for binomial", {
  set.seed(5); n <- 150; cs <- data.frame(px = runif(n), py = runif(n))
  eta <- 0.5 + sin(3 * cs$px); yp <- rpois(n, exp(eta)); yb <- rbinom(n, 1, plogis(eta - 0.5))
  hp <- quiet(cf_glm_hv(y = yp, coords = cs, family = poisson()))
  mp <- quiet(cf_glm(y = yp, coords = cs, x0 = NULL, coords0 = far, mod_hv = hp))
  expect_true(all(is.finite(mp$pred0$pred_sd)))
  hb <- quiet(cf_glm_hv(y = yb, coords = cs, family = binomial()))
  b1 <- quiet(cf_glm(y = yb, coords = cs, coords0 = far, mod_hv = hb))
  b0 <- withr::with_options(list(spcf.stage_bound = FALSE), quiet(cf_glm(y = yb, coords = cs, coords0 = far, mod_hv = hb)))
  expect_identical(b1$pred0, b0$pred0)
})

## ---- cf_dglm ----
st <- sim_spacetime(ns = 60, nt = 4, seed = 9)
hvd <- quiet(cf_dglm_hv(y = st$y, x = st$x, coords = st$coords, time = st$time))
far_st <- data.frame(px = c(0.5, 1.5, 3), py = c(0.5, 1.5, 3))

test_that("cf_dglm stage_bound keeps the point predictions and gives finite variances", {
  a <- quiet(cf_dglm(y = st$y, x = st$x, coords = st$coords, time = st$time,
                     x0 = st$x[1:3, , drop = FALSE], coords0 = far_st, time0 = rep(st$nt, 3),
                     mod_hv = hvd, se_type = "mean"))
  b <- withr::with_options(list(spcf.stage_bound = FALSE), quiet(cf_dglm(y = st$y, x = st$x, coords = st$coords, time = st$time,
                     x0 = st$x[1:3, , drop = FALSE], coords0 = far_st, time0 = rep(st$nt, 3),
                     mod_hv = hvd, se_type = "mean")))
  expect_equal(a$pred0$pred, b$pred0$pred)
  expect_true(all(is.finite(a$pred0$pred_sd) & a$pred0$pred_sd > 0))
  expect_true(is.finite(hvd$other$tau_stage) || length(hvd$other$bands) == 0)
})

test_that("cf_dglm adds the unmodeled field variance when no scale is accepted", {
  set.seed(3); ns <- 40; nt <- 3; cs <- data.frame(px = runif(ns), py = runif(ns))
  y <- rnorm(ns * nt) + rep(sin(6 * cs$px), nt)
  co <- cs[rep(seq_len(ns), nt), ]; tt <- rep(seq_len(nt), each = ns)
  h0 <- quiet(cf_dglm_hv(y = y, coords = co, time = tt)); h0$other$bands <- numeric(0)
  m1 <- quiet(cf_dglm(y = y, coords = co, time = tt, coords0 = cs[1:5, ], time0 = rep(nt, 5),
                      mod_hv = h0, se_type = "mean"))
  m0 <- withr::with_options(list(spcf.stage_bound = FALSE), quiet(cf_dglm(y = y, coords = co, time = tt, coords0 = cs[1:5, ], time0 = rep(nt, 5),
                      mod_hv = h0, se_type = "mean")))
  expect_equal(m1$pred0$pred, m0$pred0$pred)
  expect_true(all(m1$pred0$pred_sd >= m0$pred0$pred_sd))
})
