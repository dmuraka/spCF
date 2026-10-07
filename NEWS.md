# spCF 0.2.2.1

* Tests only: the `predict()` tests compared predictions made in different
  batches for exact identity. With an alternative BLAS (BLIS, OpenBLAS) the
  matrix products of different sizes are evaluated in a different order, so
  the results agree only up to rounding (about 1e-16); they are now compared
  with a tolerance of 1e-10. No change to the package code.

# spCF 0.2.2

**Effect on results:**

* `cf_lm()`, `cf_glm()`: the same accepted scales, coefficients and point
  predictions as 0.2.1. Predictive SDs and quantiles change (per-stage bound
  on the field variance, except for binomial), and so, slightly, do the
  default coefficient SEs, which use the calibrated field variance.
* `cf_dglm()`: the same accepted scales, coefficients and point predictions as
  0.2.1, unless `time0` contains time points inside the training period that
  have no observations: 0.2.1 inserted them into the AR(1) time grid, which
  changed the fit itself (in a test, by about 0.5% for predictions at the
  sample sites and 6-7% at such a time point). Predictive SDs and coefficient
  SEs change (reworked variance of the mean). `validation_MAE` in `e_summary`
  is now the mean absolute error.
* `cf_downscale()`: numerically identical to 0.2.1; the intercept row of
  `beta` is labelled `Intercept`.

## New features

* `negbin()`: negative binomial family with estimated (or fixed) dispersion
  for `cf_glm_hv()`, `cf_glm()`, `cf_dglm_hv()` and `cf_dglm()`; theta is
  re-estimated after the initial GLM and after each accepted scale.
  `poisson(link = "identity")` is also supported.
* `se_type = "prediction"` now returns a moment-matched observation predictive,
  calibrated to 95% holdout coverage, for the Gamma, inverse Gaussian,
  quasipoisson and quasibinomial families.
* `predict()` methods for `cf_lm()`, `cf_glm()` and `cf_dglm()` fits:
  prediction at new sites (and, for `cf_dglm()`, new time points) without the
  training data. The fits keep the local estimates at the
  knots of every selected scale (well under 1 MB in typical fits), so the cost
  of a prediction does not depend on the size of the training data (about 0.25
  s for 22,500 sites, for 3,000 or 100,000 training sites alike). The result
  is identical to fitting with the same sites as `coords0`.
  `predict()` returns the predictive mean, SD and quantiles at any levels
  (`probs`); without new sites it gives these at the sample sites.
* `cf_lm()`, `cf_glm()` and `cf_dglm()` gain `keep_scales` (default `TRUE`).
  With `keep_scales = FALSE` the scale-wise processes `Z`, `Z_sd`, `Z0` and
  `Z0_sd` are not kept, which makes the fitted object several times smaller;
  predictions are unchanged, and only `sp_scalewise()` needs them.

## Changes in cf_dglm() for prediction times

* Time points in `time0` that are not training time points no longer enter the
  AR(1) time grid of the fit. They are predicted from the smoothed per-knot
  states: bridged between the neighbouring training time points (the AR(1) step
  split in proportion to the time differences), or forecast / backcast by the
  time difference over the median spacing of the training time points. The fit,
  `beta_tv` and `sd_summary` therefore no longer depend on `time0` (an interior
  time point used to add a step to the AR(1) grid and so changed the fit), and
  `predict()` reproduces the predictions later. Forecasts one spacing ahead are
  unchanged.

## Smaller fitted objects

* The quantile tables `pred_q`, `pred0_q`, `pred_q_signal` and `pred0_q_signal`
  (15 columns each) are no longer stored. `mod$pred_q` and the other fields
  still return them at the same 15 levels, computed on access and identical to
  the stored tables of earlier versions.
* The prediction tables `pred` and `pred0` had character row names inherited
  from named prediction vectors, which made them about five times larger than
  their numbers; they now carry automatic row names.
* Together, a `cf_dglm()` fit of 2,000 sites x 100 time points shrinks from
  155 MB to 108 MB (41 MB with `keep_scales = FALSE`), of which 24 MB are the
  knot states kept for `predict()`, and a `cf_lm()` fit of 100,000 sites from
  90 MB to 69 MB (11 MB).

## Changes to predictive variances

* `cf_lm()` and `cf_glm()`: the calibrated variance of each spatial scale is
  now bounded by that scale's share of the field variance, so the predictive
  SD grows smoothly with the distance to the data instead of drawing rings
  around isolated sites (see Details in `?cf_lm`). Point predictions and
  coefficients are unchanged.
* `cf_dglm()`: a distance-aware field variance of the mean, an
  information-scaled calibration, and a floor on the field variance (see
  Details in `?cf_dglm`). Coverage of 95% mean intervals at held-out sites is
  close to nominal in simulations (it was 0.58-0.75). Point predictions are
  unchanged.
* The argument `sill_cap` of `cf_dglm()` is removed: the field variance is
  always capped at the marginal variance of the fitted field (as with the
  default `sill_cap = TRUE` before).
* The noise part of the `cf_lm()` coefficient standard errors (`se_method =
  "opt"`) is rescaled to a nearest-neighbour nugget estimate.

## Bug fixes

* `spCFmap()`: irregular sites with coordinates on a fine common resolution
  (e.g. integer metres, as the meuse sample sites) were taken for a lattice and
  drawn as invisible one-metre pixels; they are now filled from the nearest site.

* `spCFmap()`: irregular prediction sites are still drawn as a raster filled
  from the nearest site, but only within a circle of a common radius around
  each site, so nothing far from a prediction site is coloured. The default
  radius is 0.75 times the median distance to the nearest site (at least 1/300
  of the diagonal of the region), and a "Circle size" slider scales it. Regular
  lattices are drawn as before. For a `cf_dglm()` fit with prediction sites,
  the time slider spans the prediction time points (min(time0) to max(time0))
  instead of the whole training period, and steps by their spacing when they
  are equally spaced.

* `cf_dglm()`: `validation_MAE` in `e_summary` was the absolute mean error
  (absolute bias) instead of the mean absolute error.
* `cf_downscale()`: the intercept row of `beta` was labelled `x` or `V1`; it is
  now `Intercept`, and unnamed covariates are labelled `x1`, `x2`, ...
* `cf_lm()`: no longer warns when the holdout predictions are constant;
  `validation_R2` is then `NA`.
* The calibration of the per-scale variance no longer fails when the holdout
  moment equation has no root in its search interval.

# spCF 0.2.1

**Effect on results:** none. For the same inputs, `cf_lm()`, `cf_glm()`,
`cf_dglm()` and `cf_downscale()` reproduce 0.2.0 exactly (selection,
coefficients, SEs, predictions and calibration).

* [interface] `print()` of `cf_lm()` and `cf_glm()` fits shows tables on stdout
  again (0.2.0 sent them to stderr as deparsed vectors).
* [interface] The intercept row of `cf_glm()`'s `beta` is labelled `Intercept`
  again (it was blank in 0.2.0 with `robust_se = TRUE` or an offset).
* [interface] `cf_dglm()` keeps `x` and `x0` in `other`, so that `spCFmap()` can
  draw the covariate effect.
* [interface] `spCFmap()`: a clear error when terra/PROJ cannot build the CRS;
  keyless Esri/OSM basemaps; the covariate-effect layer is not drawn for fits
  that did not keep their covariates; irregular prediction grids are filled
  from the nearest site; colour scales of space-time layers span all time
  points.

# spCF 0.2.0

**Effect on results:**

* `cf_lm()`, `cf_glm()`: the same accepted scales, coefficients and point
  predictions as 0.1.2 (up to rounding), except binomial predicted
  probabilities, which are now calibrated (see `se_type`). Coefficient SEs,
  predictive SDs and quantiles change by default, often considerably (new
  defaults `se_method = "opt"` and `se_type = "prediction"`);
  `se_type = "mean", se_method = "classic"` reproduces 0.1.2.
* `cf_downscale()`: numerically identical to 0.1.2 (input checks and messages
  only).
* `cf_dglm()`: new.

## New features

* `cf_dglm_hv()` and `cf_dglm()`: coarse-to-fine dynamic (space-time) GLMMs,
  with a per-knot AR(1) Kalman smoother in time and optional time-varying
  coefficients (`tvc`).
* `spCFmap()`: interactive Shiny map explorer for fitted models, and a full app
  to upload data and fit models.
* `sp_scalewise()` gains `time_range` (time averages for `cf_dglm()` fits).

## Uncertainty

* [uncertainty] `cf_lm()`, `cf_glm()` (and `cf_dglm()`) gain
  `se_method = c("opt", "classic")`, default `"opt"`, for the cluster-robust
  coefficient SEs (`robust_se = TRUE`). `"classic"` is the sandwich of 0.1.2,
  which keeps the realised field in the residual; `"opt"` uses the
  field-removed noise plus the calibrated field variance, with a leave-one-out
  ceiling. Coefficients are unchanged; their SEs (above all the intercept's)
  and the coefficient part of the predictive SD change.
* [uncertainty] `se_type = c("prediction", "mean")`, default `"prediction"`:
  `pred_sd` and the quantiles describe a new observation, calibrated on the
  holdout samples of the `*_hv` fit (Gaussian: residual variance with a
  split-conformal factor; Poisson: negative-binomial predictive; binomial:
  temperature-calibrated probability). `"mean"` gives the predictive
  distribution of the mean, as in 0.1.2. The mean versions are kept as
  `pred_signal` and `pred_q_signal`; the calibration is in `other$calibration`.
* [estimates] Binomial `cf_glm()` with the default `se_type = "prediction"`:
  `pred$pred` and `pred0$pred` are the calibrated probabilities
  `plogis(qlogis(p) / T)`, with T fitted on the holdout samples.

## Computation

* [speed/memory] The local regressions of `cf_lm()` and `cf_glm()` run in a
  fused C++ kernel with a bundled k-d tree (nanoflann) instead of chunked
  neighbour lists: same results up to rounding, faster, with much lower peak
  memory.
* [selection] `cf_glm_hv()` no longer tries bandwidths below half the median
  nearest-neighbour distance of the sites. Such scales have almost no data in
  their kernels; the selection changes only if 0.1.2 had accepted one of them.
* [estimates] `cf_lm(add_learn = "lightgbm")`: the final refit uses the number
  of boosting rounds scaled to all data.

## Interface

* [interface] Input checks with informative errors in all `cf_*()` functions.
* [interface] Progress output goes through `message()`, so `suppressMessages()`
  silences it.
* [interface] `sp_scalewise()`: `bw_range` is half-open, `[min, max)`.
* [interface] `ranger` moved from Imports to Suggests; R >= 4.1.0 required.
* [interface] `cf_downscale()` keeps the coordinates in `other$coords`.

# spCF 0.1.2

**Effect on results:**

* `cf_lm()`, `cf_glm()`: coefficients and point predictions change slightly (a
  change in the local kernel weighting; about 1% for predictions and 0.1% for
  coefficients in a test), and the accepted scales can change. SEs and
  predictive SDs change substantially (new robust SEs and calibrated
  predictive variance).
* `cf_downscale()`: new.

* `cf_downscale_hv()` and `cf_downscale()`: coarse-to-fine spatial downscaling
  of areal data to fine units, with predictions that add up to the areal
  values.
* [estimates] [selection] In the local regressions of `cf_lm()` and `cf_glm()`
  (C++), the predictive variance used to weight the local estimates is
  `sigma / w` instead of `x^2 sigma / sum(w^2) + sigma / w^2`. This changes the
  combined local estimates and, through the holdout loss, the accepted scales.
* [uncertainty] `cf_lm()` and `cf_glm()` gain `robust_se = TRUE`: spatial-block
  cluster-robust coefficient SEs (also used in the coefficient part of the
  predictive SD).
* [uncertainty] The spatial part of the predictive variance is the predictive
  variance of the local estimates (instead of their coefficient variance),
  multiplied by a factor calibrated on the holdout samples (`other$tau`) and
  capped at the variance of the fitted field (not for binomial).
* [uncertainty] `cf_glm()`: `pred0_q` uses the link-scale SD (it used the
  response-scale SD on the link scale).
* [uncertainty] `cf_lm()`: predictive quantiles `pred_q` and `pred0_q`; with
  `add_learn`, they combine the model and the learner by conformalized quantile
  regression on the holdout samples.
* [interface] `e_summary` of `cf_lm()` and `cf_glm()` is computed on the holdout
  samples of the `*_hv` fit; the MAE of `cf_glm()` is the mean absolute error
  (it was the absolute mean error).
* [interface] `cf_lm_hv(add_learn = "lightgbm")`.

# spCF 0.1.1

**Effect on results:**

* `cf_lm()`: the training/validation split changes (k-means up to 30,000
  sites, seeded), and hence the accepted scales and estimates.
* `cf_glm()`: new.

* `cf_glm_hv()` and `cf_glm()`: coarse-to-fine spatial GLMMs for any `glm`
  family, with offsets and prediction sites. The scales are selected on the
  holdout deviance; after each accepted scale the GLM is refitted with the
  accumulated field as an offset, the field is centred and the intercept
  re-based; the linear predictor is clipped at +-20 for non-identity links
  (`options(spcf.l_pred_cap)`). Predictive quantiles `pred_q` and `pred0_q` are
  on the response scale.
* [selection] `cf_lm_hv()` gains `seed = 123`. The training/validation split
  uses k-means for up to 30,000 unique sites (1,000 in 0.1.0), otherwise
  random sampling, and is reproducible.
* [selection] Fine scales with more than 1,000 knots place them by seeded random
  sampling instead of k-means.
* [estimates] `cf_lm()`: a negative scaling of the fitted field is no longer
  set to zero (edge case).
* [speed/memory] The local regressions run in C++ (Rcpp), over chunks of knots.
* [interface] `cf_lm()` accepts `x = NULL` (intercept only); `Z`, `Z_sd`, `Z0`
  and `Z0_sd` are data frames with named columns.

# spCF 0.1.0

* First CRAN release: `cf_lm_hv()`, `cf_lm()` and `sp_scalewise()` for the
  Gaussian coarse-to-fine spatial model.

## Notes on this file (all versions)

Tags: [estimates] changes coefficients or point predictions for the same data
and arguments; [selection] changes the scales accepted by the `*_hv` functions
(and hence the estimates); [uncertainty] changes only SEs, predictive SDs,
intervals or quantiles; [speed/memory] computation only, same results up to
rounding; [interface] arguments, outputs or printing only. The notes on
0.1.0-0.2.1 were reconstructed from the CRAN releases.
