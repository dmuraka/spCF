## ---------------------------------------------------------------------------
## Predictive quantiles of cf_lm / cf_glm / cf_dglm fits.
##
## The fits do not store quantile tables (for large data, and above all for
## space-time data, these tables made up a large part of the fitted object);
## they keep the predictive mean and SD on the link scale together with the
## holdout calibration, from which .spcf_quantile() rebuilds the quantiles
## exactly. Users get them through mod$pred_q / mod$pred0_q (the 15 default
## levels) and, at any levels, through predict(..., probs).
## ---------------------------------------------------------------------------

.spcf_check_probs <- function(probs) {
  if (!is.numeric(probs) || !length(probs) || anyNA(probs) || any(probs <= 0 | probs >= 1))
    .spcf_stop("'probs' must be probabilities strictly between 0 and 1.")
  invisible(probs)
}

## levels of mod$pred_q (the 15 default levels)
.spcf_probs_of <- function(mod) .spcf_qs

## mod: a cf_lm / cf_glm / cf_dglm fit. sites: "sample" (rows of mod$pred) or
## "prediction" (rows of mod$pred0). type: "prediction" (the predictive
## distribution that mod$pred describes) or "signal" (that of the mean).
.spcf_quantile <- function(mod, probs = NULL, sites = c("sample", "prediction"),
                           type = c("prediction", "signal")) {
  if (!inherits(mod, c("cf_lm", "cf_glm", "cf_dglm")))
    .spcf_stop("'mod' must be a fit from cf_lm(), cf_glm() or cf_dglm().")
  sites <- match.arg(sites); type <- match.arg(type)
  if (is.null(probs)) probs <- .spcf_probs_of(mod)
  .spcf_check_probs(probs)
  o    <- .subset2(mod, "other")
  spec <- o$qspec
  if (is.null(spec)) {
    ## a fit saved by spCF <= 0.2.1 still carries its quantile tables
    nm <- paste0(if (sites == "sample") "pred_q" else "pred0_q",
                 if (type == "signal" && identical(o$se_type, "prediction")) "_signal" else "")
    Q  <- .subset2(mod, nm)
    if (is.null(Q)) return(NULL)
    j  <- match(paste0("q", probs), names(Q))
    if (anyNA(j)) .spcf_stop("This fit (spCF <= 0.2.1) stores the quantiles only at the levels ",
                             paste(sub("^q", "", names(Q)), collapse = ", "), ".")
    return(Q[, j, drop = FALSE])
  }
  if (type == "prediction" && identical(o$se_type, "prediction")) {
    q <- .spcf_obs_q(mod, sites, probs)
    if (!is.null(q)) return(q)
  }
  .spcf_signal_q(spec, sites, probs)
}

## Quantile tables of earlier versions, now rebuilt on access by .spcf_quantile()
## at the 15 default levels.
.spcf_qfields <- c(pred_q = "sample", pred0_q = "prediction",
                   pred_q_signal = "sample", pred0_q_signal = "prediction")

.spcf_get <- function(x, name) {
  v <- .subset2(x, name)
  if (!is.null(v) || !(name %in% names(.spcf_qfields)) || is.null(.subset2(x, "other")$qspec))
    return(v)
  signal <- grepl("_signal$", name)
  ## pred_q_signal existed only for fits with se_type = "prediction"
  if (signal && !identical(.subset2(x, "other")$se_type, "prediction")) return(NULL)
  .spcf_quantile(x, sites = .spcf_qfields[[name]], type = if (signal) "signal" else "prediction")
}

## Only the four quantile fields are intercepted; every other name goes to the
## default method (with its partial matching), so nothing else changes.
#' @export
`$.cf_lm`    <- function(x, name) if (name %in% names(.spcf_qfields)) .spcf_get(x, name) else NextMethod()
#' @export
`$.cf_glm`   <- function(x, name) if (name %in% names(.spcf_qfields)) .spcf_get(x, name) else NextMethod()
#' @export
`$.cf_dglm`  <- function(x, name) if (name %in% names(.spcf_qfields)) .spcf_get(x, name) else NextMethod()
#' @export
`[[.cf_lm`   <- function(x, i, ...) if (is.character(i) && length(i) == 1L && i %in% names(.spcf_qfields)) .spcf_get(x, i) else NextMethod()
#' @export
`[[.cf_glm`  <- function(x, i, ...) if (is.character(i) && length(i) == 1L && i %in% names(.spcf_qfields)) .spcf_get(x, i) else NextMethod()
#' @export
`[[.cf_dglm` <- function(x, i, ...) if (is.character(i) && length(i) == 1L && i %in% names(.spcf_qfields)) .spcf_get(x, i) else NextMethod()
