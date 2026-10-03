#' Negative binomial family with an estimated (or fixed) dispersion
#'
#' A \code{\link{family}} object for negative binomial responses, for use as
#' the \code{family} argument of \code{\link{cf_glm_hv}} and
#' \code{\link{cf_dglm_hv}}. With \code{theta = NULL} (default) the dispersion
#' parameter \eqn{\theta} (variance \eqn{\mu + \mu^2/\theta}) is estimated by
#' maximum likelihood: it is re-estimated on the training samples after the
#' initial GLM and after every accepted spatial scale, so that it reflects the
#' over-dispersion that remains once the multiscale spatial process is modeled.
#' A numeric \code{theta} keeps the dispersion fixed (equivalent to
#' \code{MASS::negative.binomial(theta)}, which is also accepted).
#'
#' @param theta \code{NULL} (estimate) or a positive number (fixed).
#' @param link Link function: \code{"log"} (default), \code{"sqrt"} or
#'   \code{"identity"}.
#'
#' @return An object of class \code{"family"}. The fitted value of
#'   \eqn{\theta} is stored in \code{family$theta} of the family returned in
#'   \code{mod_hv$other$family} (and \code{mod$other$family}).
#'
#' @seealso \code{\link{cf_glm_hv}}, \code{\link{cf_dglm_hv}}
#'
#' @examples
#' set.seed(1)
#' n      <- 300
#' coords <- cbind(runif(n), runif(n))
#' y      <- rnbinom(n, size = 2, mu = exp(1 + sin(4 * coords[, 1])))
#' mod_hv <- cf_glm_hv(y = y, coords = coords, family = negbin())
#' mod_hv$other$family$theta   # estimated dispersion
#'
#' @export
negbin <- function(theta = NULL, link = "log") {
  if (!is.null(theta) && (length(theta) != 1L || !is.finite(theta) || theta <= 0))
    .spcf_stop("'theta' must be NULL (estimate) or a single positive number.")
  if (!link %in% c("log", "sqrt", "identity"))
    .spcf_stop("negbin() supports the 'log', 'sqrt' and 'identity' links.")
  .spcf_nb_family(if (is.null(theta)) 1 else theta, link, estimate = is.null(theta))
}

## Negative binomial family object (the same form as MASS::negative.binomial,
## without the MASS dependency). `estimate` marks theta for re-estimation.
#' @keywords internal
#' @noRd
.spcf_nb_family <- function(theta, link = "log", estimate = FALSE) {
  lk <- stats::make.link(link)
  variance   <- function(mu) mu + mu^2 / theta
  validmu    <- function(mu) all(is.finite(mu)) && all(mu > 0)
  dev.resids <- function(y, mu, wt)
    2 * wt * (y * log(pmax(1, y) / mu) - (y + theta) * log((y + theta) / (mu + theta)))
  aic <- function(y, n, mu, wt, dev) {
    term <- (y + theta) * log(mu + theta) - y * log(mu) + lgamma(y + 1) -
      theta * log(theta) + lgamma(theta) - lgamma(theta + y)
    2 * sum(term * wt)
  }
  initialize <- expression({
    if (any(y < 0)) stop("negative values not allowed for the negative binomial family")
    n <- rep(1, nobs)
    mustart <- y + (y == 0) / 6
  })
  fam <- list(family = paste0("Negative Binomial(", format(round(theta, 4)), ")"),
              link = link, linkfun = lk$linkfun, linkinv = lk$linkinv,
              variance = variance, dev.resids = dev.resids, aic = aic,
              mu.eta = lk$mu.eta, initialize = initialize, validmu = validmu,
              valideta = lk$valideta, theta = theta,
              spcf_estimate_theta = isTRUE(estimate))
  fam <- structure(fam, class = "family")
  if (link %in% c("identity", "sqrt")) fam <- .spcf_positive_mean(fam)
  fam
}

#' @keywords internal
#' @noRd
.spcf_is_nb <- function(family) is.character(family$family) &&
  startsWith(family$family, "Negative Binomial")

## theta of a negative binomial family (MASS objects carry it only in the name).
#' @keywords internal
#' @noRd
.spcf_nb_theta <- function(family) {
  if (!is.null(family$theta)) return(family$theta)
  suppressWarnings(as.numeric(sub("^Negative Binomial\\((.*)\\)$", "\\1", family$family)))
}

## Maximum-likelihood theta given the current means (optionally weighted),
## bounded to [1e-3, 1e6]. Falls back to the previous value if the search fails.
#' @keywords internal
#' @noRd
.spcf_theta_ml <- function(y, mu, wt = NULL, old = NA_real_) {
  ok <- is.finite(y) & is.finite(mu) & mu > 0
  if (is.null(wt)) wt <- rep(1, length(y))
  y <- y[ok]; mu <- mu[ok]; wt <- wt[ok]
  if (length(y) < 3L) return(if (is.finite(old)) old else 1)
  ll <- function(lt) {
    th <- exp(lt)
    sum(wt * (lgamma(th + y) - lgamma(th) - lgamma(y + 1) + th * log(th) +
                y * log(mu) - (th + y) * log(th + mu)))
  }
  o  <- tryCatch(stats::optimize(ll, c(log(1e-3), log(1e6)), maximum = TRUE),
                 error = function(e) NULL)
  th <- if (is.null(o) || !is.finite(o$maximum)) old else exp(o$maximum)
  if (!is.finite(th)) th <- 1
  min(max(th, 1e-3), 1e6)
}

## Re-estimate theta of a negbin() family with theta to be estimated (other
## families, including MASS::negative.binomial with fixed theta, are returned
## unchanged). `idx` restricts the estimation to the training samples.
#' @keywords internal
#' @noRd
.spcf_nb_update <- function(family, y, mu, idx = NULL) {
  if (!isTRUE(family$spcf_estimate_theta)) return(family)
  if (!is.null(idx)) { y <- y[idx]; mu <- mu[idx] }
  th <- .spcf_theta_ml(y, mu, old = family$theta)
  .spcf_nb_family(th, family$link, estimate = TRUE)
}

## Positive-mean guard for families whose link does not keep the mean in the
## valid range (Poisson/negative binomial with the identity or sqrt link). The
## inverse link is floored at a small positive value, so the IRLS of glm() /
## glm.fit() never stops at a non-positive mean ("no valid set of coefficients
## has been found") when the spatial offset pushes the linear predictor below
## zero; the floor acts like the mu floor used for the log link.
#' @keywords internal
#' @noRd
.spcf_positive_mean <- function(family, floor = 1e-8) {
  if (isTRUE(family$spcf_positive_mean)) return(family)
  li <- family$linkinv; me <- family$mu.eta
  if (identical(family$link, "identity")) {
    family$linkinv <- function(eta) pmax(eta, floor)
    family$mu.eta  <- function(eta) ifelse(eta > floor, 1, 1e-8)
  } else if (identical(family$link, "sqrt")) {
    family$linkinv <- function(eta) pmax(eta, sqrt(floor))^2
    family$mu.eta  <- function(eta) 2 * pmax(eta, sqrt(floor))
  }
  family$valideta <- function(eta) TRUE
  family$validmu  <- function(mu) all(is.finite(mu)) && all(mu > 0)
  family$spcf_positive_mean <- TRUE
  family
}

## Prepare the user-supplied family once at the cf_*_hv entry point:
##  * accepts a family object, a family function or its name (as glm() does);
##  * poisson(identity): positive-mean guard (otherwise glm() fails to start);
##  * gaussian with a log/inverse link: positive starting values when some
##    responses are non-positive (otherwise glm() cannot find starting values).
## Other families are returned unchanged, so existing fits are not affected.
#' @keywords internal
#' @noRd
.spcf_prepare_family <- function(family) {
  if (is.character(family)) family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) family <- family()
  if (!inherits(family, "family"))
    .spcf_stop("'family' must be a family object such as poisson(), binomial() or negbin().")
  if (identical(family$family, "poisson") && identical(family$link, "identity"))
    family <- .spcf_positive_mean(family)
  if (identical(family$family, "gaussian") && family$link %in% c("log", "inverse")) {
    family$initialize <- expression({
      n <- rep.int(1, nobs)
      if (is.null(etastart) && is.null(start) && is.null(mustart)) {
        pos <- y[y > 0]
        mustart <- pmax(y, if (length(pos)) 0.1 * min(pos) else 1e-3)
      } else mustart <- y
    })
  }
  family
}
