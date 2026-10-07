// Fused per-scale operator for CF-DGLMM: neighbour-limited kernel aggregation
// -> per-knot AR(1) Kalman filter + RTS smoother -> generalized Product-of-
// Experts (gPoE) recombination, in a single pass over the sparse neighbour
// structure. Replaces the chain of Matrix sparse products + .dglm_ksmooth,
// avoiding intermediate K x T matrices and repeated sparse traversals.
//
// LAYOUT: every hot array is stored time-major (T is the leading/contiguous
// dimension) so the innermost t-loop runs over contiguous memory. The working
// panels arrive transposed as W0t, R0t (T x nL, column-major) so W0t[t + T*i]
// is contiguous in t. Aggregates/state/gPoE accumulators are T x K or T x n.
//
// Neighbour lists are CSR: for site i, its knots are idx[ptr[i] .. ptr[i+1]-1]
// with kernel weights w[...] (0-based knot indices).
#include <Rcpp.h>
#include <vector>
using namespace Rcpp;

// [[Rcpp::export]]
List dglm_scale_chunk(IntegerVector ptr, IntegerVector idx, NumericVector w,
                      NumericMatrix W0t, NumericMatrix R0t,
                      int K, double rho, double Q,
                      IntegerVector pptr, IntegerVector pidx, NumericVector pw,
                      int n0, int return_state = 0) {
  const int T  = W0t.nrow();          // time-major: rows = time
  const int nL = W0t.ncol();          //            cols = sites
  const double eps = 1e-12;
  const double *pW0 = &W0t[0], *pR0 = &R0t[0];

  // ---- 1. aggregate working residual to knots: den, Znum, Rnum (T x K) ----
  std::vector<double> den((size_t)T * K, 0.0), Znum((size_t)T * K, 0.0),
                      Rnum((size_t)T * K, 0.0);
  for (int i = 0; i < nL; ++i) {
    const double *w0i = pW0 + (size_t)T * i;
    const double *r0i = pR0 + (size_t)T * i;
    for (int nz = ptr[i]; nz < ptr[i + 1]; ++nz) {
      const int k = idx[nz];
      const double wik = w[nz], wik2 = wik * wik;
      double *dk = &den[(size_t)T * k], *zk = &Znum[(size_t)T * k],
             *rk = &Rnum[(size_t)T * k];
      for (int t = 0; t < T; ++t) {
        const double w0 = w0i[t];
        if (w0 == 0.0) continue;
        const double wW = wik * w0;
        dk[t] += wW;
        zk[t] += wW * r0i[t];
        rk[t] += wik2 * w0;
      }
    }
  }

  // ---- 2. per-knot AR(1) Kalman filter + RTS smoother -> m, P (T x K) ----
  std::vector<double> m((size_t)T * K), P((size_t)T * K);
  std::vector<double> af(T), Pf(T), ap(T), Pp(T);
  // filtered variances, kept (return_state = 1) for bridging to new time points
  std::vector<double> Pf_all(return_state ? (size_t)T * K : 0);
  const double P0 = Q / (1.0 - rho * rho);
  for (int k = 0; k < K; ++k) {
    const double *dk = &den[(size_t)T * k], *zk = &Znum[(size_t)T * k],
                 *rk = &Rnum[(size_t)T * k];
    double *mk = &m[(size_t)T * k], *Pk = &P[(size_t)T * k];
    double a = 0.0, p = P0;
    for (int t = 0; t < T; ++t) {
      const double a_pred = rho * a;
      const double p_pred = rho * rho * p + Q;
      ap[t] = a_pred; Pp[t] = p_pred;
      const double d = dk[t];
      if (d > eps) {                                  // observed knot-time -> update
        const double Zkt = zk[t] / d;
        const double Rkt = rk[t] / (d * d);
        const double Kg  = p_pred / (p_pred + Rkt);
        a = a_pred + Kg * (Zkt - a_pred);
        p = (1.0 - Kg) * p_pred;
      } else { a = a_pred; p = p_pred; }              // missing -> predict only
      af[t] = a; Pf[t] = p;
      if (return_state) Pf_all[(size_t)T * k + t] = p;
    }
    double ms = af[T - 1], ps = Pf[T - 1];
    mk[T - 1] = ms; Pk[T - 1] = ps > 1e-8 ? ps : 1e-8;
    for (int t = T - 2; t >= 0; --t) {
      double pp1 = Pp[t + 1]; if (pp1 < 1e-12) pp1 = 1e-12;
      const double G = rho * Pf[t] / pp1;
      ms = af[t] + G * (ms - ap[t + 1]);
      ps = Pf[t] + G * G * (ps - Pp[t + 1]);
      mk[t] = ms; Pk[t] = ps > 1e-8 ? ps : 1e-8;
    }
  }
  // precompute knot precision and mean*precision (T x K) for gPoE
  std::vector<double> invP((size_t)T * K), mP((size_t)T * K);
  for (size_t j = 0; j < invP.size(); ++j) { invP[j] = 1.0 / P[j]; mP[j] = m[j] * invP[j]; }

  // Variance-only path for Vd (does not affect m, P or the mean F). The knot
  // prior P0 = Q/(1-rho^2) comes from one (rho, Q) shared by all scales and is
  // generally far from the amplitude of this scale, so the posterior/prior ratio
  // P_kt/P0 carries little information. Given the scale's own prior variance
  // P0v (set below from the variance of its fitted knot means), the Kalman
  // variance recursion -- which depends only on (rho, Qv) and the knot
  // observation variances, not on the data values -- is rerun with
  // Qv = P0v (1 - rho^2). Pv (T x K) then feeds Vd.
  std::vector<double> Pv((size_t)T * K);
  std::vector<double> vPf_all(return_state ? (size_t)T * K : 0);
  auto var_path = [&](double P0v) {
    const double Qv = P0v * (1.0 - rho * rho);
    std::vector<double> vPf(T), vPp(T);
    for (int k = 0; k < K; ++k) {
      const double *dk = &den[(size_t)T * k], *rk = &Rnum[(size_t)T * k];
      double *Pk = &Pv[(size_t)T * k];
      double p = P0v;
      for (int t = 0; t < T; ++t) {
        const double p_pred = rho * rho * p + Qv; vPp[t] = p_pred;
        const double d = dk[t];
        if (d > eps) { const double Rkt = rk[t] / (d * d); p = (1.0 - p_pred / (p_pred + Rkt)) * p_pred; }
        else p = p_pred;
        vPf[t] = p;
        if (return_state) vPf_all[(size_t)T * k + t] = p;
      }
      double ps = vPf[T - 1]; Pk[T - 1] = ps > 1e-12 ? ps : 1e-12;
      for (int t = T - 2; t >= 0; --t) {
        double pp1 = vPp[t + 1]; if (pp1 < 1e-12) pp1 = 1e-12;
        const double G = rho * vPf[t] / pp1;
        ps = vPf[t] + G * G * (ps - vPp[t + 1]);
        Pk[t] = ps > 1e-12 ? ps : 1e-12;
      }
    }
  };
  // scale prior for the variance path: variance of the smoothed knot means over
  // observed knot-times (the amplitude of this scale); falls back to P0
  double P0v = P0;
  {
    double s1 = 0.0, s2 = 0.0; long nn = 0;
    for (int k = 0; k < K; ++k) for (int t = 0; t < T; ++t)
      if (den[(size_t)T * k + t] > eps) { const double v = m[(size_t)T * k + t]; s1 += v; s2 += v * v; ++nn; }
    if (nn > 1) { const double mu = s1 / nn, vv = (s2 - nn * mu * mu) / (nn - 1); if (vv > 1e-12) P0v = vv; }
  }
  var_path(P0v);

  // ---- 3. gPoE recombination over a set of sites' CSR neighbours ----
  // Generalized product of experts with kernel weights NORMALIZED to sum to one
  // (Cao & Fleet 2014): with k~_ik = w_ik / sum_k w_ik,
  //   F(i,t) = sum_k (w_ik/P_kt) m_kt / sum_k (w_ik/P_kt)   (normalizer cancels),
  //   V(i,t) = 1 / sum_k (k~_ik/P_kt) = (sum_k w_ik) / sum_k (w_ik/P_kt).
  // The mean is unchanged versus the unnormalized form; only the variance differs
  // (it no longer shrinks purely with the number of nearby knots).
  // Vd (distance-aware variance, returned alongside V; the mean F is unchanged):
  // a knot's posterior variance Pv_kt (variance path above, prior P0v) informs
  // the field at site i through the kernel correlation w_ik = exp(-d/b), so the
  // conditional variance of the site field given that knot is
  // w^2 Pv + (1 - w^2) P0v. Vd = sum_k w_ik / sum_k [w_ik / (w^2 Pv + (1-w^2) P0v)]
  // grows from Pv_kt next to a well-determined knot to P0v at the edge of the
  // kernel support, and equals P0v where no knot reaches the site. The returned
  // "P0" is P0v, so Vd / P0 is the fraction of the scale's prior that remains. V keeps the
  // original gPoE variance (fallback max(1, max V)) for backward compatibility.
  double vmx_tr = 1.0;                         // fallback V of the training recombination
  auto gpoe = [&](IntegerVector P_ptr, IntegerVector P_idx, NumericVector P_w,
                  int n, NumericMatrix &F, NumericMatrix &V, NumericMatrix &Vd,
                  double *vmx_out) {
    std::vector<double> gden((size_t)T * n, 0.0), gnum((size_t)T * n, 0.0),
                        gdd((size_t)T * n, 0.0), sumw((size_t)n, 0.0);
    for (int i = 0; i < n; ++i) {
      double *gd = &gden[(size_t)T * i], *gn = &gnum[(size_t)T * i], *g2 = &gdd[(size_t)T * i];
      double sw = 0.0;
      for (int nz = P_ptr[i]; nz < P_ptr[i + 1]; ++nz) {
        const int k = P_idx[nz];
        const double wik = P_w[nz];
        sw += wik;                              // sum of kernel weights (t-independent)
        const double *iPk = &invP[(size_t)T * k], *mPk = &mP[(size_t)T * k], *Pk = &Pv[(size_t)T * k];
        const double c2 = wik * wik > 1.0 ? 1.0 : wik * wik;
        for (int t = 0; t < T; ++t) {
          gd[t] += wik * iPk[t]; gn[t] += wik * mPk[t];
          g2[t] += wik / (c2 * Pk[t] + (1.0 - c2) * P0v);
        }
      }
      sumw[i] = sw;
    }
    double vmx = 1.0;
    for (int i = 0; i < n; ++i)
      for (int t = 0; t < T; ++t) {
        const double g = gden[(size_t)T * i + t];
        if (g > 0.0) { double v = sumw[i] / g; if (v > vmx) vmx = v; }
      }
    if (vmx_out) *vmx_out = vmx;
    for (int i = 0; i < n; ++i) {
      const double *gd = &gden[(size_t)T * i], *gn = &gnum[(size_t)T * i];
      const double sw = sumw[i];
      for (int t = 0; t < T; ++t) {
        const double g = gd[t];
        if (g > 0.0) { F(i, t) = gn[t] / g; V(i, t) = sw / g; }
        else         { F(i, t) = 0.0;       V(i, t) = vmx;    }
        const double g2 = gdd[(size_t)T * i + t];
        Vd(i, t) = g2 > 0.0 ? sw / g2 : P0v;
      }
    }
  };

  NumericMatrix Ftr(nL, T), Vtr(nL, T), Vtr_d(nL, T);
  gpoe(ptr, idx, w, nL, Ftr, Vtr, Vtr_d, &vmx_tr);
  List out = List::create(_["Ftr"] = Ftr, _["Vtr"] = Vtr, _["Vtr_d"] = Vtr_d, _["P0"] = P0v);
  if (n0 > 0) {
    NumericMatrix Fpr(n0, T), Vpr(n0, T), Vpr_d(n0, T);
    gpoe(pptr, pidx, pw, n0, Fpr, Vpr, Vpr_d, 0);
    out["Fpr"] = Fpr; out["Vpr"] = Vpr; out["Vpr_d"] = Vpr_d;
  }
  if (return_state) {
    // knot state (T x K, time-major as stored): smoothed mean m and variance P,
    // filtered variance Pf, and the variance path (smoothed Pv, filtered vPf)
    auto mk = [&](const std::vector<double> &v) {
      NumericMatrix M(T, K); std::copy(v.begin(), v.end(), M.begin()); return M; };
    out["state"] = List::create(_["m"] = mk(m), _["P"] = mk(P), _["Pf"] = mk(Pf_all),
                                _["Pv"] = mk(Pv), _["vPf"] = mk(vPf_all),
                                _["P0v"] = P0v, _["vmx"] = vmx_tr);
  }
  return out;
}

// gPoE recombination of stored knot states at sites that each carry one time
// column (tcol, 0-based, into the T' x K state matrices m, P, Pv). Repeats the
// arithmetic of the gpoe() lambda of dglm_scale_chunk for that column, so a
// site at a training time gets exactly the value of the fit.
// [[Rcpp::export]]
List dglm_gpoe_rows(IntegerVector ptr, IntegerVector idx, NumericVector w,
                    IntegerVector tcol, NumericMatrix m, NumericMatrix P,
                    NumericMatrix Pv, double P0v, double vmx) {
  const int n = tcol.size(), T = m.nrow();
  NumericVector F(n), V(n), Vd(n);
  for (int i = 0; i < n; ++i) {
    const int t = tcol[i];
    double gd = 0.0, gn = 0.0, g2 = 0.0, sw = 0.0;
    for (int nz = ptr[i]; nz < ptr[i + 1]; ++nz) {
      const int k = idx[nz];
      const double wik = w[nz];
      sw += wik;
      const size_t j = (size_t)T * k + t;
      const double iP = 1.0 / P[j], mP = m[j] * iP;
      const double c2 = wik * wik > 1.0 ? 1.0 : wik * wik;
      gd += wik * iP; gn += wik * mP;
      g2 += wik / (c2 * Pv[j] + (1.0 - c2) * P0v);
    }
    if (gd > 0.0) { F[i] = gn / gd; V[i] = sw / gd; }
    else          { F[i] = 0.0;     V[i] = vmx;     }
    Vd[i] = g2 > 0.0 ? sw / g2 : P0v;
  }
  return List::create(_["F"] = F, _["V"] = V, _["Vd"] = Vd);
}
