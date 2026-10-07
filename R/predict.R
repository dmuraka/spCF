## ---------------------------------------------------------------------------
## Prediction at new sites from a fitted model, without the training data.
##
## cf_lm() and cf_glm() keep, for every accepted scale, the state of its knots
## (position, local coefficient and its variance, local residual variance and
## the inverse sum of squared kernel weights; lwr_glm_fused_cpp with
## return_state = 1), together with the quantities of the fit that do not
## depend on the prediction sites (coefficients, scale weights, variance
## calibration). The prediction at new sites only scatters the knot states to
## those sites (lwr_scatter0_cpp) and combines them as the fit does.
## cf_lm()/cf_glm() themselves compute their prediction sites with the same
## functions, so predict() and the fitting functions agree exactly.
## ---------------------------------------------------------------------------

## Scale-wise intercept field and its predictive variance at new sites
## (columns of Z0 and Z0_pv), as lwr()/lwr_glm() return them for coords0.
.spcf_scatter_Z0 <- function(core, x0, coords0) {
  n0 <- nrow(x0); nS <- core$nS
  Z0 <- Z0_pv <- matrix(0, nrow = n0, ncol = nS)
  for (sc in core$scales) {
    st  <- sc$state
    r   <- lwr_scatter0_cpp(st$cent, st$iw, st$b, st$ibv, st$sig, st$vc,
                            as.matrix(coords0), x0, st$band, st$kernel_id, st$threshold)
    vc  <- st$vc
    b0  <- r$b_all0; pvi <- r$pv_inv_all0
    b0[, vc]  <- b0[, vc] / pvi[, vc]
    b0[, -vc] <- 0
    b0[is.nan(b0)] <- 0; b0[is.na(b0)] <- 0
    pv  <- 1 / pvi
    pv[, -vc] <- NA
    pv[is.nan(pv)] <- Inf
    b0pv <- pv[, 1]; b0pv[!is.finite(b0pv)] <- if (core$stage_bound) Inf else 0
    Z0[, sc$ii]    <- b0[, 1] - sc$zshift
    Z0_pv[, sc$ii] <- sqrt(b0pv)
  }
  list(Z0 = Z0, Z0_pv = Z0_pv)
}

## Covariates at new sites in the design of the fit (intercept + selected columns)
.spcf_design0 <- function(core, x0, n0) {
  one0 <- matrix(1, nrow = n0, ncol = 1)
  if (is.null(core$x_sel) || sum(core$x_sel) == 0) one0
  else cbind(one0, as.matrix(x0)[, core$x_sel])
}

## Calibrated field variance and scale-wise SDs at new sites
.spcf_field_var0 <- function(core, Z0_pv) {
  if (core$stage_bound) {
    Vst0 <- .spcf_stage_var(Z0_pv, core$caps, core$tau)
    list(fv = rowSums(Vst0), Z0_sd = sqrt(Vst0))
  } else {
    field_var0 <- rowSums(Z0_pv^2)
    fv0 <- pmin(core$tau * field_var0, core$sill)
    list(fv = fv0, Z0_sd = Z0_pv * sqrt(ifelse(field_var0 > 0, fv0 / field_var0, 1)))
  }
}

## Additional learner (cf_lm add_learn) at new sites: point prediction and the
## raw quantile table, as add_mod() returns them for coords0
.spcf_addlearn0 <- function(a_mod, x0, coords0) {
  a_data0 <- data.frame(x0[, -1], coords0); names(a_data0) <- a_mod$a_xname
  if (identical(a_mod$add_learn, "rf")) {
    list(pred0 = predict(a_mod$mod, data = a_data0)$predictions,
         qmat0 = predict(a_mod$mod, data = a_data0, type = "quantiles",
                         quantiles = a_mod$qlevels)$predictions)
  } else {
    X0 <- as.matrix(a_data0[, a_mod$a_xname])
    list(pred0 = predict(a_mod$pmod, X0),
         qmat0 = t(apply(sapply(a_mod$mod, predict, X0), 1, sort)))
  }
}

## cf_lm at new sites. Returns the pieces cf_lm() stores for its prediction sites.
.spcf_lm_new <- function(core, x0, coords0, a_mod = NULL) {
  n0  <- nrow(coords0)
  x0  <- .spcf_design0(core, x0, n0)
  zz  <- .spcf_scatter_Z0(core, x0, coords0)
  beta0 <- matrix(core$beta_int, nrow = n0, ncol = core$nx, byrow = TRUE)
  if (!is.null(core$w)) beta0[, 1] <- beta0[, 1] + zz$Z0 %*% core$w
  a_pred0 <- 0; qmat0 <- NULL
  if (core$a_run) { al <- .spcf_addlearn0(a_mod, x0, coords0); a_pred0 <- al$pred0; qmat0 <- al$qmat0 }
  pred0     <- rowSums(x0 * beta0) + a_pred0
  coef_var0 <- rowSums((x0 %*% core$vmat) * x0)
  fv        <- .spcf_field_var0(core, zz$Z0_pv)
  q0 <- NULL
  if (core$a_run) {
    Qtot0 <- total_qmat(pred0, sqrt(coef_var0 + fv$fv), qmat0, a_mod$qlevels, core$qlev)
    if (!is.null(core$cqr_off)) Qtot0 <- apply_cqr(Qtot0, core$qlev, core$cqr_off)
    pred0_sd <- (Qtot0[, length(core$qlev)] - Qtot0[, 1]) / (2 * stats::qnorm(core$qlev[length(core$qlev)]))
    q0 <- unname(Qtot0)
  } else {
    pred0_sd <- sqrt(coef_var0 + fv$fv)
  }
  list(x0 = x0, pred0 = as.numeric(pred0), pred0_sd = as.numeric(pred0_sd),
       lin0 = as.numeric(pred0), lin0_sd = as.numeric(pred0_sd), q0 = q0,
       Z0 = zz$Z0, Z0_pv = zz$Z0_pv, Z0_sd = fv$Z0_sd)
}

## cf_glm at new sites
.spcf_glm_new <- function(core, x0, coords0, offset0 = NULL) {
  n0  <- nrow(coords0)
  x0  <- .spcf_design0(core, x0, n0)
  if (is.null(offset0)) offset0 <- rep(0, n0)
  zz  <- .spcf_scatter_Z0(core, x0, coords0)
  family <- core$family
  l0  <- core$beta_int[1] + rowSums(zz$Z0)
  ## linear predictor of the final GLM (y ~ offset(l_pred + offset) + covariates)
  X0  <- cbind(1, x0[, -1, drop = FALSE]); cf <- core$gcoef; ok <- !is.na(cf)
  lin0 <- drop(X0[, ok, drop = FALSE] %*% cf[ok]) + (.spcf_clip_l(l0, family) + offset0)
  fv   <- .spcf_field_var0(core, zz$Z0_pv)
  lin0_sd <- sqrt(rowSums((x0 %*% core$vmat) * x0) + fv$fv)
  list(x0 = x0, offset0 = offset0, pred0 = as.numeric(family$linkinv(lin0)),
       pred0_sd = as.numeric(response_se(pred_lin = lin0, pred_lin_sd = lin0_sd, family = family)),
       lin0 = as.numeric(lin0), lin0_sd = as.numeric(lin0_sd),
       Z0 = zz$Z0, Z0_pv = zz$Z0_pv, Z0_sd = fv$Z0_sd)
}

## Observation predictive at new sites from the stored calibration (the
## out-of-sample half of .spcf_obs_predict()): calibrated mean (binomial) and SD
.spcf_obs_new <- function(o, family, mu, s_link) {
  cal <- o$calibration; fam <- o$obs_family
  H   <- .spcf_obs_helpers(family, fam); cc <- cal$scale
  switch(cal$type,
    gaussian_conformal = {
      s <- H$to_resp(s_link, mu); list(pred = NULL, sd = cc * sqrt(s^2 + cal$sigma2)) },
    poisson_negbin = {
      s <- H$to_log(s_link, mu); list(pred = NULL, sd = sqrt(mu + (mu * cc * s)^2)) },
    negbin = {
      th <- cal$theta; s <- H$to_log(s_link, mu)
      list(pred = NULL, sd = sqrt(mu + mu^2 / (1 / (1 / th + (1 + 1 / th) * (cc * s)^2)))) },
    binomial_temperature = {
      p <- stats::plogis(stats::qlogis(.spcf_clp(mu)) / cal$temperature)
      list(pred = p, sd = sqrt(p * (1 - p))) },
    {
      kind <- sub("_moment$", "", cal$type)
      Vf   <- function(m) pmax(family$variance(m), 1e-12)
      pin  <- .spcf_obs_pos(mu, kind); sr <- H$to_resp(s_link, pin)
      list(pred = NULL, sd = sqrt(.spcf_obs_varf(pin, sr, cc, kind, cal$dispersion, Vf))) })
}

## Fill the prediction-site fields of a fitted object from .spcf_lm_new() /
## .spcf_glm_new() output (used by the fitting functions and by predict()).
.spcf_put_new <- function(res, nw, coords0) {
  o <- res$other; spec <- o$qspec
  res$pred0 <- data.frame(pred = nw$pred0, pred_sd = nw$pred0_sd)
  if (isTRUE(spec$stored)) spec$q0 <- nw$q0 else { spec$lin0 <- nw$lin0; spec$lin0_sd <- nw$lin0_sd }
  o$qspec <- spec; o$n0 <- nrow(coords0); o$x0 <- nw$x0; o$coords0 <- coords0
  if (!is.null(nw$offset0)) o$offset0 <- nw$offset0
  if (!is.null(nw$time0))   o$time0   <- nw$time0
  if (isTRUE(o$keep_scales) && length(res$bands) > 0) {
    nm <- names(res$Z)
    res$Z0 <- as.data.frame(nw$Z0); res$Z0_sd <- as.data.frame(nw$Z0_sd)
    names(res$Z0) <- names(res$Z0_sd) <- nm
    if ("Z_pv" %in% names(o)) o$Z0_pv <- nw$Z0_pv
  }
  if (identical(o$se_type, "prediction") && !is.null(o$calibration)) {
    res$pred0_signal <- res$pred0
    s0 <- .spcf_signal_slink(spec, "prediction", spec$family, o$obs_family)
    ob <- .spcf_obs_new(o, spec$family, nw$pred0, s0)
    if (!is.null(ob$pred)) res$pred0$pred <- ob$pred
    res$pred0$pred_sd <- ob$sd
  }
  res$other <- o
  res
}

## ---------------------------------------------------------------------------
## predict() methods
## ---------------------------------------------------------------------------

.spcf_predict <- function(object, x0, coords0, offset0, probs, se_type, new_fun, time0 = NULL) {
  o <- .subset2(object, "other")
  fit_se <- if (is.null(o$se_type)) "mean" else o$se_type
  if (is.null(se_type)) se_type <- fit_se
  se_type <- match.arg(se_type, c("prediction", "mean"))
  if (se_type == "prediction" && fit_se != "prediction")
    .spcf_stop("The model was fitted with se_type = \"mean\"; the observation predictive (se_type = \"prediction\") is not available.")
  if (is.null(probs)) probs <- .spcf_probs_of(object)
  .spcf_check_probs(probs)
  type <- if (se_type == "prediction") "prediction" else "signal"
  pick <- function(res, sites) {
    nm <- if (sites == "sample") "pred" else "pred0"
    if (type == "signal" && fit_se == "prediction") nm <- paste0(nm, "_signal")
    out <- .subset2(res, nm)
    q   <- .spcf_quantile(res, probs = probs, sites = sites, type = type)
    if (is.null(q)) out else cbind(out, q)
  }
  if (is.null(coords0)) return(pick(object, "sample"))   # fitted values at the sample sites

  core <- o$pcore
  if (is.null(core))
    .spcf_stop("This fit has no stored knot states (fitted by spCF <= 0.2.1); refit it with the current version to use predict() at new sites.")
  has_x <- if (!is.null(core$has_x)) isTRUE(core$has_x) && length(core$x_sel) > 0
           else !is.null(core$x_sel) && sum(core$x_sel) > 0
  if (has_x && is.null(x0))
    .spcf_stop("'x0' must be provided: the model has covariates, and the prediction sites need the same ones.")
  ## the covariate columns of the fit, as a zero-row template for the checks
  xt <- if (has_x) matrix(numeric(0), 0, length(core$x_sel), dimnames = list(NULL, core$xcols)) else NULL
  if (has_x && is.null(core$xcols) && is.null(dim(x0))) xt <- numeric(0)
  .spcf_check_newdata(x = xt, x0 = if (has_x) x0 else NULL, coords0 = coords0,
                      time0 = time0, offset0 = offset0)
  nw  <- new_fun(core, x0, coords0)
  res <- .spcf_put_new(object, nw, as.matrix(coords0))
  pick(res, "prediction")
}

#' Prediction from a fitted coarse-to-fine model
#'
#' Predicts at new sites from a \code{\link{cf_lm}} or \code{\link{cf_glm}} fit.
#' The fitted object keeps, for every selected scale, the local estimates at
#' the knots of that scale; prediction only spreads these to the new sites, so
#' its cost does not depend on the size of the training data, and the training
#' data are not needed. The result is identical to fitting the model with the
#' same sites given as \code{coords0} (and \code{x0}, \code{offset0}).
#'
#' @param object A fitted model from \code{\link{cf_lm}} or \code{\link{cf_glm}}.
#' @param x0 Covariates at the prediction sites, with the same columns as
#'   \code{x} in the fit. Required when the model has covariates.
#' @param coords0 Coordinates of the prediction sites (matrix or data.frame with
#'   two columns). If \code{NULL}, the predictions at the sample sites are
#'   returned.
#' @param offset0 Offset at the prediction sites (\code{cf_glm} only; zero if
#'   \code{NULL}).
#' @param probs Probability levels of the predictive quantiles. Defaults to the
#'   levels of \code{pred_q} in the fit (0.005, 0.025, 0.05,
#'   0.1, ..., 0.9, 0.95, 0.975, 0.995).
#' @param se_type \code{"prediction"} for the predictive distribution of a new
#'   observation or \code{"mean"} for that of the mean. Defaults to the
#'   \code{se_type} of the fit; \code{"prediction"} needs a fit with
#'   \code{se_type = "prediction"}.
#' @param ... Not used.
#'
#' @return A data.frame with one row per site: the predictive mean
#'   (\code{pred}), the predictive standard deviation (\code{pred_sd}) and the
#'   predictive quantiles (\code{q<level>}, e.g. \code{q0.025}), on the response
#'   scale.
#'
#' @details With an additional learner (\code{add_learn} in
#'   \code{\link{cf_lm_hv}}), the learner's model is kept in the fit and the
#'   quantiles of the combined predictive are simulated, as in
#'   \code{\link{cf_lm}}; they then vary slightly from call to call.
#'
#' @examples
#' set.seed(1)
#' n      <- 300
#' coords <- cbind(px = runif(n), py = runif(n))
#' x      <- data.frame(x1 = rnorm(n))
#' y      <- 0.5 * x$x1 + sin(4 * coords[, 1]) + rnorm(n, sd = 0.3)
#' hv     <- cf_lm_hv(y = y, x = x, coords = coords)
#' mod    <- cf_lm(y = y, x = x, coords = coords, mod_hv = hv)
#'
#' coords0 <- cbind(px = runif(5), py = runif(5))
#' x0      <- data.frame(x1 = rnorm(5))
#' predict(mod, x0 = x0, coords0 = coords0, probs = c(0.025, 0.975))
#' @name spCF-predict
#' @rdname spCF-predict
#' @aliases predict.cf_lm predict.cf_glm predict.cf_dglm
#' @importFrom stats predict
#' @export
predict.cf_lm <- function(object, x0 = NULL, coords0 = NULL, probs = NULL,
                          se_type = NULL, ...) {
  a_mod <- .subset2(object, "other")$a_mod
  .spcf_predict(object, x0, coords0, NULL, probs, se_type,
                function(core, x0, coords0) .spcf_lm_new(core, x0, coords0, a_mod = a_mod))
}

#' @rdname spCF-predict
#' @export
predict.cf_glm <- function(object, x0 = NULL, coords0 = NULL, offset0 = NULL,
                           probs = NULL, se_type = NULL, ...) {
  .spcf_predict(object, x0, coords0, offset0, probs, se_type,
                function(core, x0, coords0) .spcf_glm_new(core, x0, coords0, offset0))
}

## ---------------------------------------------------------------------------
## cf_dglm at new sites and times
##
## The knot states are stored per training time point (smoothed mean m and
## variance P, filtered variance Pf; and the variance path Pv, vPf of the
## distance-aware variance). A new time point is
##   * a training time: the stored state, so the fit's own values result;
##   * after the last training time: the AR(1) forecast from the last state,
##     h = (t - t_last) / dt steps ahead (dt = median spacing of the training
##     times; one step reproduces the Kalman predict step of the fit);
##   * before the first training time: the stationary AR(1) backcast;
##   * between two training times: the AR(1) bridge between their smoothed
##     states, with the step split in proportion to the time differences. The
##     split uses the continuous-time (Ornstein-Uhlenbeck) embedding of the AR(1)
##     step, so the two parts compose exactly to the step of the fit, which is
##     therefore unchanged. With rho < 0 the mean is interpolated linearly
##     between whole steps (see .dglm_state_at()).
## The time-varying coefficients (random walk) are bridged and forecast likewise.
## ---------------------------------------------------------------------------

## Locate a new time t in the training times: list(kind, j, h)
.dglm_locate <- function(t, lev, dt) {
  nT <- length(lev); j <- match(t, lev)
  if (!is.na(j)) return(list(kind = "train", j = j))
  if (t > lev[nT]) return(list(kind = "after", j = nT, h = (t - lev[nT]) / dt))
  if (t < lev[1])  return(list(kind = "before", j = 1L, h = (lev[1] - t) / dt))
  j <- max(which(lev < t))
  list(kind = "between", j = j, h = (t - lev[j]) / (lev[j + 1] - lev[j]))
}

## Knot state (vectors over knots: m, P, Pv) at a new time for one scale.
## The variances depend on rho only through rho^2, so they use the continuous
## (OU) embedding with |rho| for any sign. With rho < 0 the mean alternates in
## sign from step to step and has no continuous-time version: it is then
## interpolated linearly between whole steps (the stored states for an interior
## time; the 0, 1, 2, ... step forecasts or backcasts outside the training
## period), which is continuous and exact at whole steps.
.dglm_state_at <- function(st, loc, rho, Q) {
  m <- st$m; P <- st$P; Pf <- st$Pf; Pv <- st$Pv; vPf <- st$vPf; P0v <- st$P0v
  P0 <- Q / (1 - rho * rho); j <- loc$j; ar <- abs(rho); neg <- rho < 0
  ## mean coefficient of an h-step forecast (linear between whole steps if rho < 0)
  rmean <- function(h) if (!neg) rho^h else {
    fl <- floor(h); f <- h - fl; (1 - f) * rho^fl + f * rho^(fl + 1) }
  out <- switch(loc$kind,
    train = list(m = m[j, ], P = P[j, ], Pv = Pv[j, ]),
    after = {
      r <- ar^loc$h
      list(m = rmean(loc$h) * m[j, ], P = r * r * Pf[j, ] + Q * ((1 - r * r) / (1 - rho * rho)),
           Pv = r * r * vPf[j, ] + P0v * (1 - r * r)) },
    before = {
      r <- ar^loc$h
      list(m = rmean(loc$h) * m[j, ], P = P0 + r * r * (P[j, ] - P0), Pv = P0v + r * r * (Pv[j, ] - P0v)) },
    between = {
      h  <- loc$h; r1 <- ar^h; r2 <- ar^(1 - h)
      bridge <- function(Pf0, Ps1, Qx) {
        Pps <- r1 * r1 * Pf0 + Qx * (1 - r1 * r1) / (1 - rho * rho)
        Ppn <- pmax(r2 * r2 * Pps + Qx * (1 - r2 * r2) / (1 - rho * rho), 1e-12)
        Gs  <- r2 * Pps / Ppn
        list(Pps = Pps, Gs = Gs, P = Pps + Gs * Gs * (Ps1 - Ppn))
      }
      mp <- bridge(Pf[j, ], P[j + 1, ], Q)
      vp <- bridge(vPf[j, ], Pv[j + 1, ], P0v * (1 - rho * rho))
      mm <- if (neg) (1 - h) * m[j, ] + h * m[j + 1, ] else {
        ## filtered mean at j, recovered from the RTS step of the fit
        Pp1 <- pmax(rho * rho * Pf[j, ] + Q, 1e-12); G <- rho * Pf[j, ] / Pp1
        den <- 1 - G * rho
        af  <- ifelse(den > 1e-10, (m[j, ] - G * m[j + 1, ]) / pmax(den, 1e-10), m[j, ])
        r1 * af + mp$Gs * (m[j + 1, ] - r2 * r1 * af) }
      list(m = mm, P = mp$P, Pv = vp$P) })
  out$P  <- pmax(out$P, 1e-8)
  out$Pv <- pmax(out$Pv, 1e-12)
  out
}

## Time-varying coefficients (random walk, drift diag(q)) at a new time:
## list(beta = d-vector, V = d x d)
.dglm_tv_at <- function(tv, loc) {
  j <- loc$j; d <- ncol(tv$beta); Qm <- diag(rep_len(tv$q, d), d)
  switch(loc$kind,
    train  = list(beta = tv$beta[j, ], V = tv$V[[j]]),
    after  = list(beta = tv$beta[j, ], V = tv$V[[j]] + loc$h * Qm),
    before = list(beta = tv$beta[j, ], V = tv$V[[j]] + loc$h * Qm),
    between = {
      Pps <- tv$Pf[[j]] + loc$h * Qm; Ppn <- tv$Pf[[j]] + Qm
      G   <- Pps %*% solve(Ppn)
      list(beta = drop(tv$af[[j]] + G %*% (tv$beta[j + 1, ] - tv$af[[j]])),
           V = Pps + G %*% (tv$V[[j + 1]] - Ppn) %*% t(G)) })
}

## Neighbour lists (CSR) of the knots of one scale for each row of 'coords',
## computed once per distinct location and expanded to the rows
.dglm_nbr_rows <- function(coords, knots, band, kernel) {
  key <- paste(coords[, 1], coords[, 2], sep = "\r"); u <- !duplicated(key)
  lk  <- match(key, key[u])
  nb  <- .dglm_nbr(coords[u, , drop = FALSE], knots, band, kernel)
  len <- diff(nb$ptr)[lk]; st <- nb$ptr[lk]
  sel <- unlist(lapply(seq_along(lk), function(i) if (len[i]) st[i] + seq_len(len[i]) else integer(0)),
                use.names = FALSE)
  list(ptr = as.integer(c(0L, cumsum(len))), idx = nb$idx[sel], w = nb$w[sel])
}

.spcf_dglm_new <- function(core, x0, coords0, time0, offset0 = NULL) {
  coords0 <- as.matrix(coords0); n0 <- nrow(coords0)
  if (is.null(offset0)) offset0 <- rep(0, n0)
  X0 <- if (!isTRUE(core$has_x)) matrix(1, n0, 1)
        else cbind(1, as.matrix(x0)[, core$x_sel, drop = FALSE])
  lev <- core$lev; dt <- if (length(lev) >= 2) stats::median(diff(lev)) else 1
  tu  <- sort(unique(time0)); tk0 <- match(time0, tu)
  locs <- lapply(tu, .dglm_locate, lev = lev, dt = dt)
  family <- core$family; nS <- core$nS
  Z0 <- Z0_sd <- matrix(0, n0, max(nS, 1L)); Rd0 <- matrix(1, n0, max(nS, 1L))
  for (k in seq_len(nS)) {
    ks <- core$kstates[[k]]; st <- ks$state
    sts <- lapply(locs, .dglm_state_at, st = st, rho = core$rho, Q = core$Q)
    M  <- do.call(rbind, lapply(sts, `[[`, "m"));  P <- do.call(rbind, lapply(sts, `[[`, "P"))
    Pv <- do.call(rbind, lapply(sts, `[[`, "Pv"))
    nb <- .dglm_nbr_rows(coords0, ks$knots, ks$band, core$kernel)
    g  <- dglm_gpoe_rows(nb$ptr, nb$idx, nb$w, as.integer(tk0 - 1L),
                         matrix(M, nrow = length(tu)), matrix(P, nrow = length(tu)),
                         matrix(Pv, nrow = length(tu)), st$P0v, st$vmx)
    Z0[, k] <- g$F - ks$zmean; Z0_sd[, k] <- sqrt(g$V); Rd0[, k] <- g$Vd / st$P0v
  }
  f0 <- if (nS > 0) rowSums(Z0) else rep(0, n0)
  ## time-varying coefficients at the new times
  has_tv <- length(core$tv_cols) > 0 && !is.null(core$tv)
  tvpart0 <- tvvar0 <- rep(0, n0)
  if (has_tv) {
    tvs <- lapply(locs, .dglm_tv_at, tv = core$tv)
    X0tv <- X0[, core$tv_cols, drop = FALSE]
    B  <- do.call(rbind, lapply(tvs, `[[`, "beta"))
    tvpart0 <- rowSums(X0tv * B[tk0, , drop = FALSE])
    tvvar0  <- pmax(vapply(seq_len(n0), function(i)
      drop(X0tv[i, ] %*% tvs[[tk0[i]]]$V %*% X0tv[i, ]), numeric(1)), 0)
  }
  ## final GLM (y ~ offset(field + tv part + offset) + constant covariates)
  ncv <- length(core$const_cov)
  Xg0 <- if (ncv > 0) cbind(1, X0[, core$const_cov, drop = FALSE]) else matrix(1, n0, 1)
  off0 <- .dglm_clip_l(f0, family) + tvpart0 + offset0
  cf <- core$gcoef; ok <- !is.na(cf)
  lin0 <- drop(Xg0[, ok, drop = FALSE] %*% cf[ok]) + off0
  fv0 <- if (core$stage_on) rowSums(.spcf_stage_var_ratio(Rd0, core$caps, core$tau_stage))
         else if (core$v_unmod > 0) rep(core$v_unmod, n0)
         else pmin(core$tau * rowSums(Z0_sd^2), core$sill)
  lin0_sd <- sqrt(pmax(rowSums((Xg0 %*% core$vmat) * Xg0) + tvvar0 + fv0, 0))
  list(x0 = x0, offset0 = offset0, time0 = time0,
       pred0 = as.numeric(family$linkinv(lin0)),
       pred0_sd = as.numeric(abs(family$mu.eta(lin0)) * lin0_sd),
       lin0 = as.numeric(lin0), lin0_sd = as.numeric(lin0_sd),
       Z0 = Z0[, seq_len(nS), drop = FALSE], Z0_sd = Z0_sd[, seq_len(nS), drop = FALSE])
}

#' @rdname spCF-predict
#' @param time0 Time points of the prediction sites (\code{cf_dglm} only), one
#'   per row of \code{coords0}. They may be training time points, time points
#'   between them (bridged between the smoothed states of the neighbouring
#'   training times), or time points before or after the training period
#'   (AR(1) backcast or forecast, the number of steps being the time difference
#'   over the median spacing of the training times).
#' @export
predict.cf_dglm <- function(object, x0 = NULL, coords0 = NULL, time0 = NULL, offset0 = NULL,
                            probs = NULL, se_type = NULL, ...) {
  if (!is.null(coords0) && is.null(time0))
    .spcf_stop("'time0' must be supplied together with 'coords0': every prediction site needs a time point.")
  .spcf_predict(object, x0, coords0, offset0, probs, se_type,
                function(core, x0, coords0) .spcf_dglm_new(core, x0, coords0, time0, offset0),
                time0 = time0)
}
