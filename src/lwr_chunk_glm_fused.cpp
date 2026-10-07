// Fused serial kernel for lwr_glm: builds a nanoflann kd-tree over `coords`
// (and over `coords0` when prediction sites are given) ONCE, then for each
// knot does (radius search -> local GLM -> scatter-add) inline, WITHOUT
// materialising the neighbour lists in R. Numerically identical to the
// frNN + lwr_chunk_glm_cpp path (radius search finds the same neighbours;
// tiny ~1e-13 differences only from sqrt(squared-dist) vs direct distance).
// Much lower peak memory (no ~O(N*avg_nb) neighbour list) and a bit faster.
// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include "nanoflann.h"
#include <vector>
#include <cmath>
#include <cstdint>
using namespace Rcpp;

struct PC2 {
  const double* pts; std::size_t N;
  inline std::size_t kdtree_get_point_count() const { return N; }
  inline double kdtree_get_pt(std::size_t i, std::size_t d) const { return pts[i + d * N]; }
  template<class BBOX> bool kdtree_get_bbox(BBOX&) const { return false; }
};
typedef nanoflann::KDTreeSingleIndexAdaptor<
          nanoflann::L2_Simple_Adaptor<double, PC2>, PC2, 2> KDTree2;

static inline double kfun_f(double d, double band, int kid) {
  if (kid == 2) return std::exp(-(d * d) / (band * band));
  return std::exp(-d / band);
}

// [[Rcpp::export]]
List lwr_glm_fused_cpp(
    NumericMatrix coords,        // n x 2  (training points; tree)
    NumericMatrix coords_cent,   // n_knot x 2 (knot query locations)
    NumericVector resid,         // n
    NumericVector w_obs,         // n
    NumericMatrix x,             // n x nx
    IntegerVector id_train,      // n (0/1)
    NumericMatrix B_var,         // n_knot x nx
    IntegerVector vc_cols,       // 1-based columns with a varying coefficient
    double band, int kernel_id, double threshold, int is_lm,
    SEXP coords0_sexp,           // n0 x 2 or NULL
    SEXP x0_sexp,                // n0 x nx or NULL
    int return_state = 0) {      // 1: also return the knot state used by lwr_scatter0_cpp

  const int n      = x.nrow();
  const int nx     = x.ncol();
  const int n_knot = coords_cent.nrow();
  const int n_vc   = vc_cols.size();
  const bool has0  = !Rf_isNull(coords0_sexp);
  const double rad2 = threshold * threshold;

  const double* cp   = coords.begin();
  const double* cc   = coords_cent.begin();
  const double* xp   = x.begin();
  const double* rp   = resid.begin();
  const double* wp   = w_obs.begin();
  const int*    idt  = id_train.begin();
  const double* Bv   = B_var.begin();
  const int*    vcp  = vc_cols.begin();

  PC2 cloud{cp, (std::size_t)n};
  KDTree2 tree(2, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10));
  tree.buildIndex();

  int n0 = 0;
  NumericMatrix x0_mat;
  const double* x0p = 0;
  std::vector<double> c0buf;          // owns coords0 for the tree's lifetime
  PC2 cloud0{0, 0};
  KDTree2* tree0 = 0;
  if (has0) {
    NumericMatrix c0(coords0_sexp);
    x0_mat = NumericMatrix(x0_sexp);
    n0 = c0.nrow();
    c0buf.assign(c0.begin(), c0.end());
    cloud0.pts = c0buf.data(); cloud0.N = (std::size_t)n0;
    tree0 = new KDTree2(2, cloud0, nanoflann::KDTreeSingleIndexAdaptorParams(10));
    tree0->buildIndex();
    x0p = x0_mat.begin();
  }

  NumericMatrix b_all(n, nx), bv_inv_all(n, nx), pv_inv_all(n, nx), b_old(n_knot, nx);
  NumericMatrix b_all0, bv_inv_all0, pv_inv_all0;
  if (has0) { b_all0 = NumericMatrix(n0, nx); bv_inv_all0 = NumericMatrix(n0, nx); pv_inv_all0 = NumericMatrix(n0, nx); }
  double* ba = b_all.begin(); double* bvv = bv_inv_all.begin(); double* pvv = pv_inv_all.begin(); double* bo = b_old.begin();
  double* ba0 = has0 ? b_all0.begin() : 0; double* bv0 = has0 ? bv_inv_all0.begin() : 0; double* pv0 = has0 ? pv_inv_all0.begin() : 0;

  // Knot state (return_state = 1): everything the scatter to a prediction site
  // needs, so that predictions at new sites can be made later without the
  // training data. NaN marks a knot (st_iw) or a knot/column (st_sig) that
  // contributed nothing.
  NumericVector st_iw;  NumericMatrix st_b, st_ibv, st_sig;
  if (return_state) {
    st_iw  = NumericVector(n_knot, NA_REAL);
    st_b   = NumericMatrix(n_knot, n_vc); st_ibv = NumericMatrix(n_knot, n_vc);
    st_sig = NumericMatrix(n_knot, n_vc);
    std::fill(st_sig.begin(), st_sig.end(), NA_REAL);
  }

  std::vector<nanoflann::ResultItem<uint32_t, double> > mt, mt0;
  std::vector<double> wei, wei0;
  nanoflann::SearchParameters sprm; sprm.sorted = false;

  for (int k = 0; k < n_knot; ++k) {
    double q[2] = { cc[k], cc[k + (std::size_t)n_knot] };
    mt.clear();
    const std::size_t m = tree.radiusSearch(q, rad2, mt, sprm);
    if (m == 0) continue;
    if (wei.size() < m) wei.resize(m);

    int m_hv = 0; double wxy = 0.0, wxxw = 0.0;
    for (std::size_t i = 0; i < m; ++i) {
      const double wker = kfun_f(std::sqrt(mt[i].second), band, kernel_id);
      wei[i] = wker;
      const int sidx = (int)mt[i].first;
      if (idt[sidx]) {
        const double ww = wker * wker, w_o = wp[sidx];
        wxy += ww * w_o * rp[sidx]; wxxw += ww * w_o; ++m_hv;
      }
    }
    if (m_hv <= 5 || wxxw <= 0.0) continue;
    const double b_sel0 = wxy / wxxw;
    for (int j = 0; j < nx; ++j) bo[k + (std::size_t)j * n_knot] = b_sel0;

    std::size_t m0 = 0;
    if (has0) {
      mt0.clear();
      m0 = tree0->radiusSearch(q, rad2, mt0, sprm);
      if (m0 > 0) {
        if (wei0.size() < m0) wei0.resize(m0);
        for (std::size_t i = 0; i < m0; ++i)
          wei0[i] = kfun_f(std::sqrt(mt0[i].second), band, kernel_id);
      }
    }

    for (int vi = 0; vi < n_vc; ++vi) {
      const int j = vcp[vi] - 1;
      const double* xj = xp + (std::size_t)j * n;
      double sigma = 0.0;
      if (is_lm) {
        // LM kernel: variance over ALL neighbours (train+test), unweighted,
        // divided by (m - 1). (m > m_hv > 5, so m - 1 >= 5.)
        for (std::size_t i = 0; i < m; ++i) {
          const int sidx = (int)mt[i].first;
          const double rs = rp[sidx] - xj[sidx] * b_sel0;
          const double v  = wei[i] * rs;
          sigma += v * v;
        }
        sigma /= (m - 1);
      } else {
        // GLM kernel: variance over TRAIN neighbours only, IRLS-weighted,
        // divided by (m_hv - 1).
        for (std::size_t i = 0; i < m; ++i) {
          const int sidx = (int)mt[i].first;
          if (!idt[sidx]) continue;
          const double rs = rp[sidx] - xj[sidx] * b_sel0;
          const double v  = wei[i] * rs;
          sigma += wp[sidx] * v * v;
        }
        if (m_hv <= 1) continue;
        sigma /= (m_hv - 1);
      }
      const double lambda    = sigma / Bv[k + (std::size_t)j * n_knot];
      const double wxxw_lam  = wxxw + lambda;
      const double b_sel_val = wxy / wxxw_lam;
      const double bv_sel    = sigma / wxxw_lam;
      const double inv_bv    = 1.0 / bv_sel;
      const double inv_wxxw  = 1.0 / wxxw;
      if (return_state) {
        st_iw[k] = inv_wxxw;
        st_b(k, vi) = b_sel_val; st_ibv(k, vi) = inv_bv; st_sig(k, vi) = sigma;
      }
      double* baj = ba + (std::size_t)j * n; double* bvj = bvv + (std::size_t)j * n; double* pvj = pvv + (std::size_t)j * n;
      for (std::size_t i = 0; i < m; ++i) {
        const int sidx = (int)mt[i].first;
        const double wk = wei[i], ws = wk * wk, xv = xj[sidx];
        const double pv_sel = (xv * xv * inv_wxxw) * sigma + sigma / wk;
        const double wei2   = ws / pv_sel;
        baj[sidx] += wei2 * b_sel_val; bvj[sidx] += wei2 * inv_bv; pvj[sidx] += wei2;
      }
      if (has0 && m0 > 0) {
        const double* x0j = x0p + (std::size_t)j * n0;
        double* ba0j = ba0 + (std::size_t)j * n0; double* bv0j = bv0 + (std::size_t)j * n0; double* pv0j = pv0 + (std::size_t)j * n0;
        for (std::size_t i = 0; i < m0; ++i) {
          const int sidx0 = (int)mt0[i].first;
          const double w0 = wei0[i], w0s = w0 * w0, xv0 = x0j[sidx0];
          const double pv_sel0 = (xv0 * xv0 * inv_wxxw) * sigma + sigma / w0;
          const double wei2_0  = w0s / pv_sel0;
          ba0j[sidx0] += wei2_0 * b_sel_val; bv0j[sidx0] += wei2_0 * inv_bv; pv0j[sidx0] += wei2_0;
        }
      }
    }
  }
  if (tree0) delete tree0;

  List out = List::create(_["b_all"]=b_all,_["bv_inv_all"]=bv_inv_all,_["pv_inv_all"]=pv_inv_all,_["b_old"]=b_old);
  if (has0) {
    out["b_all0"] = b_all0; out["bv_inv_all0"] = bv_inv_all0; out["pv_inv_all0"] = pv_inv_all0;
  }
  if (return_state) {
    out["state"] = List::create(_["iw"]=st_iw, _["b"]=st_b, _["ibv"]=st_ibv, _["sig"]=st_sig);
  }
  return out;
}

// Scatter of a stored knot state (lwr_glm_fused_cpp with return_state = 1) to
// prediction sites. Repeats, knot by knot and in the same order, the has0 block
// of lwr_glm_fused_cpp, so the accumulators are identical to those that the
// fused kernel returns when it is given the same prediction sites.
// [[Rcpp::export]]
List lwr_scatter0_cpp(
    NumericMatrix coords_cent,   // n_knot x 2
    NumericVector iw,            // n_knot: 1 / sum of the knot's squared kernel weights (NaN: inactive)
    NumericMatrix b,             // n_knot x n_vc: local coefficient
    NumericMatrix ibv,           // n_knot x n_vc: 1 / its variance
    NumericMatrix sig,           // n_knot x n_vc: local residual variance (NaN: inactive)
    IntegerVector vc_cols,       // 1-based columns with a varying coefficient
    NumericMatrix coords0,       // n0 x 2
    NumericMatrix x0,            // n0 x nx
    double band, int kernel_id, double threshold) {

  const int n_knot = coords_cent.nrow();
  const int n_vc   = vc_cols.size();
  const int n0     = coords0.nrow();
  const int nx     = x0.ncol();
  const double rad2 = threshold * threshold;
  const double* cc  = coords_cent.begin();
  const double* x0p = x0.begin();
  const int*    vcp = vc_cols.begin();

  NumericMatrix b_all0(n0, nx), bv_inv_all0(n0, nx), pv_inv_all0(n0, nx);
  double* ba0 = b_all0.begin(); double* bv0 = bv_inv_all0.begin(); double* pv0 = pv_inv_all0.begin();
  if (n0 == 0 || n_knot == 0)
    return List::create(_["b_all0"]=b_all0,_["bv_inv_all0"]=bv_inv_all0,_["pv_inv_all0"]=pv_inv_all0);

  std::vector<double> c0buf(coords0.begin(), coords0.end());
  PC2 cloud0{c0buf.data(), (std::size_t)n0};
  KDTree2 tree0(2, cloud0, nanoflann::KDTreeSingleIndexAdaptorParams(10));
  tree0.buildIndex();

  std::vector<nanoflann::ResultItem<uint32_t, double> > mt0;
  std::vector<double> wei0;
  nanoflann::SearchParameters sprm; sprm.sorted = false;

  for (int k = 0; k < n_knot; ++k) {
    const double inv_wxxw = iw[k];
    if (ISNAN(inv_wxxw)) continue;
    double q[2] = { cc[k], cc[k + (std::size_t)n_knot] };
    mt0.clear();
    const std::size_t m0 = tree0.radiusSearch(q, rad2, mt0, sprm);
    if (m0 == 0) continue;
    if (wei0.size() < m0) wei0.resize(m0);
    for (std::size_t i = 0; i < m0; ++i)
      wei0[i] = kfun_f(std::sqrt(mt0[i].second), band, kernel_id);
    for (int vi = 0; vi < n_vc; ++vi) {
      const double sigma = sig(k, vi);
      if (ISNAN(sigma)) continue;
      const int j = vcp[vi] - 1;
      const double b_sel_val = b(k, vi), inv_bv = ibv(k, vi);
      const double* x0j = x0p + (std::size_t)j * n0;
      double* ba0j = ba0 + (std::size_t)j * n0; double* bv0j = bv0 + (std::size_t)j * n0; double* pv0j = pv0 + (std::size_t)j * n0;
      for (std::size_t i = 0; i < m0; ++i) {
        const int sidx0 = (int)mt0[i].first;
        const double w0 = wei0[i], w0s = w0 * w0, xv0 = x0j[sidx0];
        const double pv_sel0 = (xv0 * xv0 * inv_wxxw) * sigma + sigma / w0;
        const double wei2_0  = w0s / pv_sel0;
        ba0j[sidx0] += wei2_0 * b_sel_val; bv0j[sidx0] += wei2_0 * inv_bv; pv0j[sidx0] += wei2_0;
      }
    }
  }
  return List::create(_["b_all0"]=b_all0,_["bv_inv_all0"]=bv_inv_all0,_["pv_inv_all0"]=pv_inv_all0);
}
