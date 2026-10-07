## ---- Per-stage bounded spatial-process predictive variance (cf_lm, cf_glm) ----
## The field predictive variance is a sum over the committed stages of each
## stage's pv-based predictive variance pv_r (eq.(10)). Summed as is, it has two
## defects at sparse sites: a stage whose kernel support no knot with data
## reaches contributes nothing (pv_r is non-finite and was set to 0), and near the
## edge of that support pv_r is inflated by the tiny weights of distant knots. On
## a map this alternates saturated and empty rings around isolated sites.
## Rule (default, stage_bound = TRUE): the calibrated variance of stage r is
##   min(tau * pv_r, kappa * s_r^2),   kappa = var(sum_r Z_r) / sum_r s_r^2,
## where s_r^2 is the variance of stage r's fitted field over the sample sites,
## pv_r = Inf outside its support, and tau solves
##   mean_val sum_r min(tau * pv_r, kappa s_r^2) = holdout MSE - noise.
## Near the data a stage contributes its calibrated variance; where the data say
## nothing about it, its share of the marginal field variance (the stage caps sum
## to the sill). kappa is needed because the stages are strongly positively
## correlated (sum_r s_r^2 is typically a tenth of var(sum_r Z_r)); tau is applied
## before the bound because it is set near the data (usually < 1) and would
## otherwise shrink the level far from it as well.
.spcf_stage_caps <- function(Z, sill){
  caps <- apply(Z, 2, function(z){ v <- stats::var(z); if(is.finite(v)) v else 0 })
  if(is.finite(sill) && sum(caps) > 0) caps <- caps * sill / sum(caps)
  caps
}
.spcf_stage_var <- function(P, caps, tau){
  if(is.null(P)) return(NULL)
  if(!ncol(P)) return(matrix(0, nrow(P), 0L))
  pmin(tau * P^2, matrix(caps, nrow(P), length(caps), byrow = TRUE))
}
## Solve the (weighted) moment equation for tau_raw; 1e-6 when the
## noise-removed error is not positive, 100 when it exceeds the sill.
.spcf_stage_tau_raw <- function(P, caps, num, w = rep(1, nrow(P))){
  g   <- function(t) sum(w * rowSums(.spcf_stage_var(P, caps, t))) / sum(w)
  top <- sum(caps)
  if(!ncol(P) || top <= 0) return(1)
  if(!is.finite(num) || num <= 0) return(1e-6)
  if(num >= top) return(100)
  if(g(1e-10) >= num) return(1e-10)   # no sign change at the lower end
  exp(stats::uniroot(function(lt) g(exp(lt)) - num, c(log(1e-10), log(1e10)), tol = 1e-8)$root)
}

## ---- cf_dglm: stage variance from a posterior/prior ratio (information scaling) ----
## For cf_dglm the per-scale variance is supplied as r = Vd / P0v in [0, 1] (the
## fraction of the scale's prior variance left after the data; dglm_chunk.cpp).
## The holdout factor tau scales the information the data carry rather than the
## variance: with prior variance c_r (the stage cap, caps summing to the sill) and
## data precision (1/r - 1)/c_r, scaling the data precision by 1/tau gives
##   var_r = c_r * tau r / (1 + r (tau - 1)),
## which is ~ c_r tau r near the data and equals c_r where the data say nothing
## (r = 1), whatever tau. (A variance multiplier min(tau r, 1) c_r would leave the
## level far from the data at tau c_r, below the sill when tau < 1.)
.spcf_stage_var_ratio <- function(R, caps, tau){
  if(is.null(R)) return(NULL)
  if(!ncol(R)) return(matrix(0, nrow(R), 0L))
  r <- pmin(pmax(R, 0), 1)
  sweep(tau * r / (1 + r * (tau - 1)), 2, caps, `*`)
}
.spcf_stage_tau_raw_ratio <- function(R, caps, num, w = rep(1, nrow(R))){
  g   <- function(t) sum(w * rowSums(.spcf_stage_var_ratio(R, caps, t))) / sum(w)
  top <- g(1e10)
  if(!ncol(R) || top <= 0) return(1)
  if(!is.finite(num) || num <= 0) return(1e-6)
  if(num >= top) return(100)
  if(g(1e-10) >= num) return(1e-10)   # no sign change at the lower end
  exp(stats::uniroot(function(lt) g(exp(lt)) - num, c(log(1e-10), log(1e10)), tol = 1e-8)$root)
}
