## Per-stage bounded predictive variance (stage_bound, internal_stage_var.R).

d  <- sim_spatial(n = 150, seed = 11)
hv <- quiet(cf_lm_hv(y = d$y, x = d$x, coords = d$coords))
## prediction sites from inside the data cloud to far outside it
far <- data.frame(px = c(0.5, 0.5, 1.2, 2, 4), py = c(0.5, 0.5, 1.2, 2, 4))
x0  <- data.frame(v1 = 0, v2 = 0.5)[rep(1, nrow(far)), ]

test_that("stage_bound leaves point predictions and coefficients unchanged", {
  a <- quiet(cf_lm(y = d$y, x = d$x, coords = d$coords, x0 = x0, coords0 = far, mod_hv = hv))
  b <- quiet(cf_lm(y = d$y, x = d$x, coords = d$coords, x0 = x0, coords0 = far, mod_hv = hv,
                   stage_bound = FALSE))
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
  b0 <- quiet(cf_glm(y = yb, coords = cs, coords0 = far, mod_hv = hb, stage_bound = FALSE))
  expect_identical(b1$pred0, b0$pred0)
})
