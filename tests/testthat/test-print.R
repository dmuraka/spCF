## print() of every fitted object must show tables on stdout. print.cf_lm once
## used message(format(df)), which flattened each table into deparsed character
## vectors (c("...", ...)); check every class and the edge cases that change the
## shape of the tables (no covariates, no accepted scale, count families).

expect_tables <- function(obj, must = character(0)) {
  err <- NULL
  out <- withCallingHandlers(capture.output(print(obj)),
                             message = function(m) { err <<- c(err, conditionMessage(m)); invokeRestart("muffleMessage") })
  expect_null(err)                                  # nothing printed to stderr
  expect_false(any(grepl('c("', out, fixed = TRUE)))
  expect_true(any(grepl("^Call:", out)))
  for (s in must) expect_true(any(grepl(s, out, fixed = TRUE)), info = s)
  for (k in intersect(c("beta", "sd_summary", "e_summary"), names(obj)))
    expect_true(is.data.frame(obj[[k]]) || is.matrix(obj[[k]]), info = k)
  invisible(out)
}

test_that("cf_lm and cf_lm_hv print tables, also without covariates or scales", {
  d  <- sim_spatial(n = 120)
  hv <- quiet(cf_lm_hv(y = d$y, x = d$x$v1, coords = d$coords))
  expect_tables(hv)
  expect_tables(quiet(cf_lm(y = d$y, x = d$x$v1, coords = d$coords, mod_hv = hv)),
                c("Intercept", "validation_R2", "residuals"))
  set.seed(5); yn <- rnorm(120)
  hv <- quiet(cf_lm_hv(y = yn, coords = d$coords))
  m  <- quiet(cf_lm(y = yn, coords = d$coords, mod_hv = hv))
  expect_tables(m, c("Intercept", "validation_RMSE"))
})

test_that("cf_glm prints tables for every family", {
  d   <- sim_spatial(n = 120)
  eta <- 0.3 * d$x$v1 + 0.6 * d$field
  set.seed(6)
  fams <- list(list(poisson(),  rpois(120, exp(eta))),
               list(binomial(), rbinom(120, 1, plogis(eta))),
               list(Gamma("log"), rgamma(120, 2, 2 / exp(eta))),
               list(negbin(),   rnbinom(120, size = 2, mu = exp(eta))))
  for (f in fams) {
    hv <- quiet(cf_glm_hv(y = f[[2]], x = d$x, coords = d$coords, family = f[[1]]))
    expect_tables(hv)
    expect_tables(quiet(cf_glm(y = f[[2]], x = d$x, coords = d$coords, mod_hv = hv)),
                  c("Intercept", "validation_Pseudo-R2"))
  }
})

test_that("cf_dglm and cf_downscale print tables, also without covariates", {
  s  <- sim_spacetime(ns = 40, nt = 4)
  hv <- quiet(cf_dglm_hv(y = s$y, coords = s$coords, time = s$time))
  expect_tables(hv)
  expect_tables(quiet(cf_dglm(y = s$y, coords = s$coords, time = s$time, mod_hv = hv)),
                c("Intercept", "validation_MAE"))
  a  <- sim_areal()
  hv <- quiet(cf_downscale_hv(Y = a$Y, coords = a$coords, agg_id = a$agg_id))
  expect_tables(hv)
  expect_tables(quiet(cf_downscale(Y = a$Y, coords = a$coords, agg_id = a$agg_id, mod_hv = hv)),
                c("Intercept", "validation_R2"))
})

test_that("cf_dglm validation_MAE is the mean absolute error", {
  ## it was abs(mean(error)), the absolute bias, which can be far below RMSE
  s  <- sim_spacetime(ns = 40, nt = 4)
  hv <- quiet(cf_dglm_hv(y = s$y, x = s$x, coords = s$coords, time = s$time))
  m  <- quiet(cf_dglm(y = s$y, x = s$x, coords = s$coords, time = s$time, mod_hv = hv))
  vt <- setdiff(seq_along(s$y), hv$id_train)
  e  <- s$y[vt] - m$pred$pred[vt]
  expect_equal(m$e_summary$value[3], mean(abs(e)))
})

test_that("cf_downscale labels the intercept and unnamed covariates", {
  a  <- sim_areal()
  hv <- quiet(cf_downscale_hv(Y = a$Y, x = as.matrix(unname(a$x)), coords = a$coords, agg_id = a$agg_id))
  m  <- quiet(cf_downscale(Y = a$Y, x = as.matrix(unname(a$x)), coords = a$coords, agg_id = a$agg_id, mod_hv = hv))
  expect_equal(rownames(m$beta), c("Intercept", "x1"))
  hv <- quiet(cf_downscale_hv(Y = a$Y, x = a$x, coords = a$coords, agg_id = a$agg_id))
  m  <- quiet(cf_downscale(Y = a$Y, x = a$x, coords = a$coords, agg_id = a$agg_id, mod_hv = hv))
  expect_equal(rownames(m$beta), c("Intercept", "v1"))
})

test_that("cf_lm reports NA validation_R2 without a warning when predictions are constant", {
  d  <- sim_spatial(n = 120)
  set.seed(5); yn <- rnorm(120)
  hv <- quiet(cf_lm_hv(y = yn, coords = d$coords))
  if (length(hv$other$bands) == 0L) {
    expect_no_warning(m <- suppressMessages(cf_lm(y = yn, coords = d$coords, mod_hv = hv)))
    expect_true(is.na(m$e_summary$value[1]))
  } else succeed()
})
