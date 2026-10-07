## ---------------------------------------------------------------------------
## Observation (data-distribution) predictive with holdout calibration.
##
## Converts the SIGNAL predictive (mean uncertainty, as carried in pred_q) into
## the OBSERVATION predictive for a new data point, and calibrates it on the
## cf_*_hv holdout (out-of-fold) samples so that it is verifiable on data:
##   gaussian : N(mean, signal_var + sigma^2), split-conformal SD scaling
##   poisson  : NegBin count predictive (Poisson x lognormal-lambda),
##              holdout scaling of the mean-uncertainty component
##   binomial : Bernoulli(p_calibrated); p from holdout temperature scaling,
##              pred_sd = sqrt(p(1-p)) (interval coverage is degenerate for
##              binary, so calibration is on the probability)
##
## Works purely from (pred, pred_q) + the family + the holdout info in mod_hv,
## so it is shared unchanged by cf_lm / cf_glm / cf_dglm and does not touch
## their internal variance construction. Returns replacement pred_sd / pred_q
## (and, for binomial, a calibrated pred); the caller keeps the signal versions
## in separate fields so nothing is lost (non-destructive).
## ---------------------------------------------------------------------------
.spcf_obs_predict <- function(family, y, mod_hv,
                              pred_in, s_in,
                              pred_out = NULL, s_out = NULL) {
  fam <- .spcf_obs_fam(family, y)
  clp <- .spcf_clp
  H   <- .spcf_obs_helpers(family, fam)
  lk  <- H$lk; to_resp <- H$to_resp; to_log <- H$to_log
  ## s_in / s_out: link-scale SD of the signal, recovered from the signal
  ## quantiles by .spcf_slink() (the quantiles themselves are not stored)
  if (is.null(s_in)) return(NULL)                    # nothing to build from

  ## out-of-fold holdout predictions + observed y (from cf_*_hv)
  idt     <- mod_hv$id_train
  hvp     <- mod_hv$other$pred
  val     <- if (!is.null(idt)) setdiff(seq_along(y), idt) else integer(0)
  haveval <- length(val) >= 10 && !is.null(hvp) && all(is.finite(hvp[val]))

  out <- list(type = fam)

  if (fam == "gaussian") {
    s_in  <- to_resp(s_in, pred_in)
    if (!is.null(s_out)) s_out <- to_resp(s_out, pred_out)
    sig2 <- if (haveval) mean((y[val] - hvp[val])^2) else mean((y - pred_in)^2)
    obs_in  <- sqrt(s_in^2 + sig2)
    obs_out <- if (!is.null(s_out)) sqrt(s_out^2 + sig2) else NULL
    cc <- 1
    if (haveval) {
      sv <- sqrt(s_in[val]^2 + sig2)
      cc <- as.numeric(stats::quantile(abs(y[val] - hvp[val]) / sv, 0.95,
                                       names = FALSE)) / 1.96
      if (!is.finite(cc) || cc <= 0) cc <- 1
    }
    out$pred_sd <- cc * obs_in
    if (!is.null(obs_out)) out$pred0_sd <- cc * obs_out
    out$calib <- list(type = "gaussian_conformal", sigma2 = sig2, scale = cc)

  } else if (fam == "poisson") {
    s_in  <- to_log(s_in, pred_in)
    if (!is.null(s_out)) s_out <- to_log(s_out, pred_out)
    cc <- 1
    if (haveval) {
      sv <- s_in[val]; lamv <- pmax(hvp[val], 1e-8)
      cov_at <- function(c) {
        sz <- 1 / pmax((c * sv)^2, 1e-8)
        mean(y[val] >= stats::qnbinom(.025, size = sz, mu = lamv) &
             y[val] <= stats::qnbinom(.975, size = sz, mu = lamv))
      }
      grid <- seq(0.05, 1.5, 0.05)
      cc <- grid[which.min(abs(vapply(grid, cov_at, numeric(1)) - 0.95))]
    }
    nb_sd <- function(mu, s) sqrt(mu + (mu * cc * s)^2)
    out$pred_sd <- nb_sd(pred_in, s_in)
    if (!is.null(s_out)) out$pred0_sd <- nb_sd(pred_out, s_out)
    out$calib <- list(type = "poisson_negbin", scale = cc)

  } else if (fam == "negbin") {
    ## negative binomial observation noise (theta) mixed with the lognormal
    ## uncertainty of the mean: Var = mu + mu^2/theta + (1+1/theta) mu^2 s^2,
    ## i.e. a negative binomial with 1/size = 1/theta + (1+1/theta) s^2.
    ## The mean-uncertainty scale is holdout-calibrated as for Poisson.
    th    <- .spcf_nb_theta(family)
    s_in  <- to_log(s_in, pred_in)
    if (!is.null(s_out)) s_out <- to_log(s_out, pred_out)
    size_of <- function(c, s) 1 / (1 / th + (1 + 1 / th) * (c * s)^2)
    cc <- 1
    if (haveval) {
      sv <- s_in[val]; lamv <- pmax(hvp[val], 1e-8)
      cov_at <- function(c) {
        sz <- size_of(c, sv)
        mean(y[val] >= stats::qnbinom(.025, size = sz, mu = lamv) &
             y[val] <= stats::qnbinom(.975, size = sz, mu = lamv))
      }
      grid <- seq(0.05, 1.5, 0.05)
      cc <- grid[which.min(abs(vapply(grid, cov_at, numeric(1)) - 0.95))]
    }
    nb_sd <- function(mu, s) sqrt(mu + mu^2 / size_of(cc, s))
    out$pred_sd <- nb_sd(pred_in, s_in)
    if (!is.null(s_out)) out$pred0_sd <- nb_sd(pred_out, s_out)
    out$calib <- list(type = "negbin", theta = th, scale = cc)

  } else if (fam == "binomial") {
    Tt <- 1
    if (haveval) {
      mv <- stats::qlogis(clp(hvp[val]))
      Tt <- tryCatch(stats::optimize(function(T) {
        pc <- clp(stats::plogis(mv / T))
        -mean(y[val] * log(pc) + (1 - y[val]) * log(1 - pc))
      }, c(0.3, 3))$minimum, error = function(e) 1)
    }
    pcal <- function(mu) stats::plogis(stats::qlogis(clp(mu)) / Tt)
    p_in <- pcal(pred_in)
    out$pred <- p_in; out$pred_sd <- sqrt(p_in * (1 - p_in))
    if (!is.null(pred_out)) {
      p_out <- pcal(pred_out); out$pred0 <- p_out; out$pred0_sd <- sqrt(p_out * (1 - p_out))
    }
    out$calib <- list(type = "binomial_temperature", temperature = Tt)
    out$binary <- TRUE                               # interval coverage degenerate
  } else {
    ## Other families (Gamma, inverse.gaussian, quasipoisson, quasibinomial with
    ## proportions, and any other glm family): moment-matched observation
    ## predictive. The observation noise is phi * V(mu) (family variance
    ## function, dispersion phi) and the mean uncertainty is the response-scale
    ## signal SD s. phi is the holdout Pearson moment after removing the signal
    ## variance (in-sample Pearson when there is no holdout); the scale c of the
    ## mean uncertainty is chosen, as for Poisson, so that the holdout 95%
    ## interval coverage is closest to 0.95.
    ##   Gamma            : Gamma with Var = phi (mu^2 + m) + m,        m = (c mu s_log)^2
    ##   inverse.gaussian : inverse Gaussian with Var = phi (mu^3 + 3 mu m) + m
    ##   quasipoisson     : negative binomial (Poisson if Var <= mu), Var = phi mu + m
    ##   quasibinomial    : Beta for proportions, Var = phi p(1-p) + m
    ##   other            : normal, Var = phi V(mu) + m
    kind <- .spcf_obs_kind(fam)
    Vf  <- function(mu) pmax(family$variance(mu), 1e-12)
    pos <- function(mu) .spcf_obs_pos(mu, kind)
    pin  <- pos(pred_in); pout <- if (!is.null(pred_out)) pos(pred_out) else NULL
    sr_in  <- to_resp(s_in, pin)
    sr_out <- if (!is.null(s_out)) to_resp(s_out, pout) else NULL
    phi <- if (haveval) {
      hv <- pos(hvp[val]); mean(((y[val] - hv)^2 - sr_in[val]^2) / Vf(hv))
    } else mean((y - pin)^2 / Vf(pin))
    if (!is.finite(phi) || phi <= 0) phi <- 1e-8
    varf <- function(mu, sr, c) .spcf_obs_varf(mu, sr, c, kind, phi, Vf)
    qf <- function(a, mu, v) .spcf_obs_qf(a, mu, v, kind)
    cc <- 1
    if (haveval) {
      hv <- pos(hvp[val])
      cov_at <- function(c) {
        v <- varf(hv, sr_in[val], c)
        mean(y[val] >= qf(.025, hv, v) & y[val] <= qf(.975, hv, v))
      }
      grid <- seq(0, 1.5, 0.05)
      cc <- grid[which.min(abs(vapply(grid, cov_at, numeric(1)) - 0.95))]
    }
    v_in <- varf(pin, sr_in, cc)
    out$pred_sd <- sqrt(v_in)
    if (!is.null(sr_out)) out$pred0_sd <- sqrt(varf(pout, sr_out, cc))
    out$calib <- list(type = paste0(kind, "_moment"), dispersion = phi, scale = cc)
  }
  out
}

## Inverse Gaussian distribution function and quantiles (mean mu, shape lam),
## used by the observation predictive. Quantiles by vectorized bisection on the
## log scale, so the cost is linear in the number of points.
.spcf_pinvgauss <- function(x, mu, lam) {
  a <- sqrt(lam / x)
  stats::pnorm(a * (x / mu - 1)) +
    exp(2 * lam / mu + stats::pnorm(-a * (x / mu + 1), log.p = TRUE))
}
.spcf_qinvgauss <- function(p, mu, lam, iter = 80L) {
  sdl <- sqrt(log1p(mu / lam))                     # log-scale SD of a moment-matched lognormal
  lo <- log(mu) - 12 * sdl - 1; hi <- log(mu) + 12 * sdl + 1
  for (it in seq_len(iter)) {
    md <- (lo + hi) / 2
    f  <- .spcf_pinvgauss(exp(md), mu, lam)
    up <- !is.na(f) & f < p
    lo[up] <- md[up]; hi[!up] <- md[!up]
  }
  exp((lo + hi) / 2)
}

## Apply .spcf_obs_predict output onto a cf_* result's prediction fields,
## preserving the signal versions as *_signal (non-destructive). Quantiles are
## not stored: .spcf_quantile() rebuilds them from other$qspec and the calibration.
.spcf_apply_obs <- function(res, ob) {
  if (is.null(ob)) return(res)
  res$pred_signal   <- res$pred
  if (!is.null(res$pred0))   res$pred0_signal   <- res$pred0
  ## exact-name access: ob$pred0 would partially match ob$pred0_sd
  ## in-sample
  if (!is.null(ob[["pred"]]))    res$pred$pred       <- ob[["pred"]]   # binomial: calibrated prob
  if (!is.null(ob[["pred_sd"]])) res$pred$pred_sd    <- ob[["pred_sd"]]
  ## out-of-sample
  if (!is.null(res$pred0)) {
    if (!is.null(ob[["pred0"]]))    res$pred0$pred    <- ob[["pred0"]]
    if (!is.null(ob[["pred0_sd"]])) res$pred0$pred_sd <- ob[["pred0_sd"]]
  }
  res$other$se_type      <- "prediction"
  res$other$calibration  <- ob$calib
  res$other$obs_family   <- ob$type
  res$other$binary_pred  <- isTRUE(ob$binary)
  res
}

## ---------------------------------------------------------------------------
## Helpers shared by .spcf_obs_predict() (fit time) and .spcf_quantile() (on
## demand), so that quantiles rebuilt later are identical to the ones the fit
## used to store.
## ---------------------------------------------------------------------------
.spcf_qs  <- c(0.005, 0.025, 0.05, seq(0.1, 0.9, 0.1), 0.95, 0.975, 0.995)
.spcf_clp <- function(p) pmin(pmax(p, 1e-8), 1 - 1e-8)

.spcf_obs_fam <- function(family, y) {
  fam <- if (.spcf_is_nb(family)) "negbin" else family$family
  ## quasibinomial with 0/1 responses is calibrated like binomial (temperature)
  if (identical(fam, "quasibinomial") && all(y %in% c(0, 1))) fam <- "binomial"
  fam
}

## link function used to recover the link-scale signal SD from the quantiles,
## and the delta-method conversions of a link-scale SD to the response scale and
## to the log scale of the mean (exact identities for the identity and log links)
.spcf_obs_helpers <- function(family, fam) {
  lnk <- family$link; clp <- .spcf_clp
  lk  <- switch(paste(fam, lnk),
                "gaussian identity" = function(p) p,
                "poisson log"       = function(p) log(pmax(p, 1e-8)),
                "binomial logit"    = function(p) stats::qlogis(clp(p)),
                function(p) family$linkfun(
                  if (fam %in% c("binomial", "quasibinomial")) clp(p)
                  else if (lnk %in% c("log", "sqrt", "inverse", "1/mu^2")) pmax(p, 1e-8)
                  else p))
  list(lk = lk,
       to_resp = function(s, mu) if (is.null(s) || identical(lnk, "identity")) s else
         s * abs(family$mu.eta(lk(mu))),
       to_log  = function(s, mu) if (is.null(s) || identical(lnk, "log")) s else
         s * abs(family$mu.eta(lk(mu))) / pmax(mu, 1e-8))
}

## link-scale signal SD recovered from the signal 2.5% and 97.5% quantiles
.spcf_slink <- function(q025, q975, lk) (lk(q975) - lk(q025)) / (2 * 1.96)

.spcf_obs_kind <- function(fam) switch(fam, Gamma = "gamma", inverse.gaussian = "invgauss",
                                       quasipoisson = "count", quasibinomial = "beta", "normal")
.spcf_obs_pos  <- function(mu, kind) if (kind %in% c("gamma", "invgauss", "count")) pmax(mu, 1e-8) else
  if (kind == "beta") .spcf_clp(mu) else mu
.spcf_obs_varf <- function(mu, sr, c, kind, phi, Vf) {
  m <- (c * sr)^2
  switch(kind, gamma = phi * (mu^2 + m) + m, invgauss = phi * (mu^3 + 3 * mu * m) + m,
         phi * Vf(mu) + m)
}
.spcf_obs_qf <- function(a, mu, v, kind) switch(kind,
  gamma    = stats::qgamma(a, shape = mu^2 / v, rate = mu / v),
  invgauss = .spcf_qinvgauss(a, mu, mu^3 / v),
  count    = { z <- stats::qpois(a, mu); i <- v > mu
               z[i] <- stats::qnbinom(a, size = mu[i]^2 / (v[i] - mu[i]), mu = mu[i]); z },
  beta     = { k <- pmax(mu * (1 - mu) / pmin(v, mu * (1 - mu) * (1 - 1e-6)) - 1, 1e-6)
               stats::qbeta(a, mu * k, (1 - mu) * k) },
  stats::qnorm(a, mu, sqrt(v)))

## quantile table (one column per level) from a per-level function
.spcf_qtab <- function(probs, n, f) {
  d <- as.data.frame(matrix(vapply(probs, f, numeric(n)), nrow = n))
  names(d) <- paste0("q", probs); d
}

## Signal quantiles from the stored predictive parameters (qspec): Gaussian on
## the link scale, or the stored table for cf_lm with an additional learner.
.spcf_signal_q <- function(spec, sites, probs) {
  if (isTRUE(spec$stored)) {
    Q <- if (sites == "sample") spec$q else spec$q0
    if (is.null(Q)) return(NULL)
    j <- match(signif(probs, 10), signif(spec$probs, 10))
    if (anyNA(j))
      .spcf_stop("With an additional learner (add_learn) the quantiles are available only at the levels ",
                 paste(spec$probs, collapse = ", "), ".")
    d <- as.data.frame(Q[, j, drop = FALSE]); names(d) <- paste0("q", probs); return(d)
  }
  lin <- if (sites == "sample") spec$lin else spec$lin0
  sd  <- if (sites == "sample") spec$lin_sd else spec$lin0_sd
  if (is.null(lin)) return(NULL)
  d <- as.data.frame(spec$family$linkinv(lin + outer(sd, stats::qnorm(probs), "*")))
  names(d) <- paste0("q", probs); d
}

## link-scale signal SD at the sample or prediction sites, as .spcf_obs_predict()
## expects it
.spcf_signal_slink <- function(spec, sites, family, fam) {
  q <- .spcf_signal_q(spec, sites, c(0.025, 0.975))
  if (is.null(q)) return(NULL)
  .spcf_slink(q[["q0.025"]], q[["q0.975"]], .spcf_obs_helpers(family, fam)$lk)
}

## Observation-predictive quantiles rebuilt from the calibration stored by
## .spcf_apply_obs() (NULL for binary responses, whose intervals are degenerate)
.spcf_obs_q <- function(mod, sites, probs) {
  o    <- .subset2(mod, "other"); cal <- o$calibration; fam <- o$obs_family
  spec <- o$qspec; family <- spec$family
  if (is.null(cal) || is.null(fam) || isTRUE(o$binary_pred)) return(NULL)
  ps <- .subset2(mod, if (sites == "sample") "pred_signal" else "pred0_signal")
  if (is.null(ps)) return(NULL)
  mu <- ps$pred; n <- length(mu)
  s  <- .spcf_signal_slink(spec, sites, family, fam)
  H  <- .spcf_obs_helpers(family, fam); cc <- cal$scale
  switch(cal$type,
    gaussian_conformal = {
      s <- H$to_resp(s, mu); obs <- sqrt(s^2 + cal$sigma2)
      d <- as.data.frame(mu + outer(cc * obs, stats::qnorm(probs)))
      names(d) <- paste0("q", probs); d },
    poisson_negbin = {
      s <- H$to_log(s, mu); sz <- 1 / pmax((cc * s)^2, 1e-8)
      .spcf_qtab(probs, n, function(a) stats::qnbinom(a, size = sz, mu = mu)) },
    negbin = {
      th <- cal$theta; s <- H$to_log(s, mu); sz <- 1 / (1 / th + (1 + 1 / th) * (cc * s)^2)
      .spcf_qtab(probs, n, function(a) stats::qnbinom(a, size = sz, mu = mu)) },
    {
      kind <- sub("_moment$", "", cal$type)
      Vf   <- function(m) pmax(family$variance(m), 1e-12)
      pin  <- .spcf_obs_pos(mu, kind); sr <- H$to_resp(s, pin)
      v    <- .spcf_obs_varf(pin, sr, cc, kind, cal$dispersion, Vf)
      .spcf_qtab(probs, n, function(a) .spcf_obs_qf(a, pin, v, kind)) })
}
