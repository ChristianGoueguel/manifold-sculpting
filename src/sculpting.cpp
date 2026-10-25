// Manifold Sculpting (Gashler, Ventura and Martinez, NIPS 2007).
//
// The algorithm follows the paper and the reference implementation of its
// authors (GManifoldSculpting in the Waffles library, which is in the public
// domain):
//
//  1. The k nearest neighbours of each point are found, and for each point p,
//     neighbour n and the neighbour m of n that is the most colinear with p,
//     the distance |p - n| and the angle (p, n, m) are recorded.
//  2. The data are rotated onto their principal components. The leading
//     `ndim` components are the preserved dimensions, the others are the
//     scaled dimensions.
//  3. Each iteration multiplies the scaled dimensions by the squishing rate,
//     scales the preserved dimensions up until the average neighbour distance
//     is restored, then visits the points breadth-first from a random seed and
//     moves each one, in the preserved dimensions only, by hill climbing on
//     the error between the current and the original distances and angles.
//  4. Once the scaled dimensions are negligible, they are dropped.
//
// The scaled dimensions are only ever multiplied by a common factor, so their
// contribution to every squared distance and dot product is that factor
// squared times its value at the start. These contributions are computed once
// and the scaled dimensions are never stored, which is how all the dimensions
// are kept (as the paper asks) at the cost of the preserved ones.

#include <RcppArmadillo.h>

#include <algorithm>
#include <cmath>
#include <deque>
#include <limits>
#include <utility>
#include <vector>

namespace {

const double kPi = 3.14159265358979323846;

// The relationship between a point p and one of its neighbours n.
struct Link {
  int nbr;        // index of n
  int slot;       // position of m in the neighbour list of n (-1: none)
  double dist0;   // original distance |p - n|
  double theta0;  // original angle (p, n, m), divided by pi
  double jsq;     // |p - n|^2 in the scaled dimensions, at scale 1
  double jdot;    // (p - n).(m - n) in the scaled dimensions, at scale 1
};

// The angle between two vectors, divided by pi, from their dot product and
// squared norms. Degenerate (zero-length) vectors give a right angle.
inline double angle01(double dot, double sq_a, double sq_b) {
  const double denom = sq_a * sq_b;
  if (dot == 0.0 || !(denom > 0.0)) {
    return 0.5;
  }
  const double c = std::max(-1.0, std::min(1.0, dot / std::sqrt(denom)));
  return std::acos(c) / kPi;
}

// 0-based, row-major copy of a 1-based n x k index matrix.
std::vector<int> flat_index(const Rcpp::IntegerMatrix& nn) {
  const int n = nn.nrow();
  const int k = nn.ncol();
  std::vector<int> out(static_cast<size_t>(n) * k);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < k; ++j) {
      out[static_cast<size_t>(i) * k + j] = nn(i, j) - 1;
    }
  }
  return out;
}

std::vector<double> flat_distance(const Rcpp::NumericMatrix& nd) {
  const int n = nd.nrow();
  const int k = nd.ncol();
  std::vector<double> out(static_cast<size_t>(n) * k);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < k; ++j) {
      out[static_cast<size_t>(i) * k + j] = nd(i, j);
    }
  }
  return out;
}

// For a point q with neighbour v (a training point), finds the neighbour of v
// that forms the angle (q, v, m) closest to pi. Returns its slot in the
// neighbour list of v, or -1, and sets `theta` to the angle divided by pi.
// `skip` is a point that cannot be m (q itself, for a training point).
int most_colinear(const double* q, int v, double sq_qv, const arma::mat& xt,
                  const std::vector<int>& idx, const std::vector<double>& dist,
                  int k, int skip, double& theta) {
  const int p = static_cast<int>(xt.n_rows);
  const double* xv = xt.colptr(v);
  int best_slot = -1;
  theta = 0.0;
  double best = -1.0;
  for (int s = 0; s < k; ++s) {
    const size_t vs = static_cast<size_t>(v) * k + s;
    const int b = idx[vs];
    if (b == skip) continue;
    const double* xb = xt.colptr(b);
    double dot = 0.0;
    for (int t = 0; t < p; ++t) {
      dot += (q[t] - xv[t]) * (xb[t] - xv[t]);
    }
    const double a = angle01(dot, sq_qv, dist[vs] * dist[vs]);
    if (a > best) {
      best = a;
      best_slot = s;
    }
  }
  if (best_slot >= 0) theta = best;
  return best_slot;
}

// In one dimension, all the angles of an unfolded manifold are pi. Waffles
// sets the original angles to pi in this case, which this follows.
inline double target_angle(double theta, int ndim) {
  if (ndim < 2) return theta > 0.0 ? 1.0 : 0.0;
  return theta;
}

class Sculptor {
 public:
  Sculptor(int n, int d, int k, std::vector<double> y, std::vector<Link> links,
           double squishing_rate)
      : n_(n), d_(d), k_(k), y_(std::move(y)), links_(std::move(links)),
        cycle_(n, -1), rate_(squishing_rate), scale_(1.0), scale2_(1.0), pass_(0) {
    double sum = 0.0;
    for (const Link& l : links_) sum += l.dist0;
    ave_dist_ = sum / static_cast<double>(links_.size());
    inv_ave_ = 1.0 / std::max(ave_dist_, 1e-300);
    eta_ = ave_dist_;

    // Bring the starting configuration to the scale of the data (a no-op for
    // the principal components, which preserve the distances).
    const double avg = average_distance();
    if (avg > 0.0) {
      const double f = ave_dist_ / avg;
      for (double& v : y_) v *= f;
      for (Link& l : links_) {
        l.jsq *= f * f;
        l.jdot *= f * f;
      }
    }
  }

  double scale() const { return scale_; }
  double learning_rate() const { return eta_; }
  double ave_dist() const { return ave_dist_; }
  const std::vector<double>& y() const { return y_; }

  // One iteration; returns the total error of the adjusted points.
  double pass() {
    // Step 4a: squish the scaled dimensions, then scale the preserved ones up
    // until the average neighbour distance is restored.
    if (scale_ > 0.0) {
      scale_ *= rate_;
      if (scale_ <= 0.001) scale_ = 0.0;
      scale2_ = scale_ * scale_;

      double avg = average_distance();
      for (int guard = 0; avg < ave_dist_ && guard < 1000; ++guard) {
        const double up = 1.0 / rate_;
        for (double& v : y_) v *= up;
        const double next = average_distance();
        if (!(next > avg)) break;  // degenerate preserved dimensions
        avg = next;
      }
    }

    // Step 4b: restore the relationships, breadth-first from a random seed
    // point, which is left in place. Points that the traversal does not reach
    // (points that are nobody's neighbour, or other components of the
    // neighbourhood graph) are visited afterwards, so that every point is
    // adjusted at each iteration.
    long steps = 0;
    double total = 0.0;
    const int seed = random_index(n_);
    traverse(seed, false, steps, total);
    const int offset = random_index(n_);
    for (int t = 0; t < n_; ++t) {
      const int i = (offset + t) % n_;
      if (cycle_[i] != pass_) traverse(i, true, steps, total);
    }

    // Learning rate: about one productive hill-climbing sweep per point.
    if (steps < n_) {
      eta_ *= 0.87;
    } else {
      eta_ /= 0.91;
    }

    if (pass_ % 20 == 0) center();
    ++pass_;
    return total;
  }

  void center() {
    for (int t = 0; t < d_; ++t) {
      double mean = 0.0;
      for (int i = 0; i < n_; ++i) mean += y_[static_cast<size_t>(i) * d_ + t];
      mean /= n_;
      for (int i = 0; i < n_; ++i) y_[static_cast<size_t>(i) * d_ + t] -= mean;
    }
  }

 private:
  static int random_index(int n) {
    const int i = static_cast<int>(R::unif_rand() * n);
    return std::min(std::max(i, 0), n - 1);
  }

  const double* row(int i) const { return &y_[static_cast<size_t>(i) * d_]; }

  double average_distance() const {
    double sum = 0.0;
    for (int i = 0; i < n_; ++i) {
      const double* yi = row(i);
      for (int j = 0; j < k_; ++j) {
        const Link& l = links_[static_cast<size_t>(i) * k_ + j];
        const double* yv = row(l.nbr);
        double sq = 0.0;
        for (int t = 0; t < d_; ++t) {
          const double a = yi[t] - yv[t];
          sq += a * a;
        }
        sum += std::sqrt(sq + scale2_ * l.jsq);
      }
    }
    return sum / static_cast<double>(links_.size());
  }

  // Equation 1 of the paper, as weighted in Waffles: the distance term is
  // normalised by the average neighbour distance, the angle term only
  // penalises angles that have closed, and the terms of the neighbours that
  // have not been adjusted yet in this iteration are down-weighted.
  double point_error(int i) const {
    const double* yi = row(i);
    double err = 0.0;
    for (int j = 0; j < k_; ++j) {
      const Link& l = links_[static_cast<size_t>(i) * k_ + j];
      const int v = l.nbr;
      const double* yv = row(v);

      double sq_a = 0.0;
      for (int t = 0; t < d_; ++t) {
        const double a = yi[t] - yv[t];
        sq_a += a * a;
      }
      sq_a += scale2_ * l.jsq;
      double dd = (std::sqrt(sq_a) - l.dist0) * inv_ave_;

      double dth = 0.0;
      if (l.slot >= 0) {
        const Link& lm = links_[static_cast<size_t>(v) * k_ + l.slot];
        const double* ym = row(lm.nbr);
        double dot = 0.0;
        double sq_b = 0.0;
        for (int t = 0; t < d_; ++t) {
          const double a = yi[t] - yv[t];
          const double b = ym[t] - yv[t];
          dot += a * b;
          sq_b += b * b;
        }
        dot += scale2_ * l.jdot;
        sq_b += scale2_ * lm.jsq;
        dth = std::min(0.0, angle01(dot, sq_a, sq_b) - l.theta0);
      }

      if (cycle_[v] != pass_) {
        dd *= 0.4;
        dth *= 0.4;
      }
      err += dd * dd + dth * dth;
    }
    return err;
  }

  // Hill climbing in the preserved dimensions. Returns the number of sweeps
  // that improved the error and sets `err` to the final error.
  int adjust_point(int i, double& err) {
    double* yi = &y_[static_cast<size_t>(i) * d_];
    double base = point_error(i);
    // The random step size helps the points get away from each other.
    const double step = eta_ * (0.6 + 0.4 * R::unif_rand());
    int sweeps = 0;
    bool progress = true;
    for (; progress && sweeps < 30; ++sweeps) {
      progress = false;
      for (int t = 0; t < d_; ++t) {
        yi[t] += step;
        double e = point_error(i);
        if (e >= base) {
          yi[t] -= 2.0 * step;
          e = point_error(i);
        }
        if (e >= base) {
          yi[t] += step;
        } else {
          base = e;
          progress = true;
        }
      }
    }
    err = base;
    return sweeps - 1;
  }

  void traverse(int start, bool adjust_start, long& steps, double& total) {
    queue_.clear();
    queue_.push_back(start);
    while (!queue_.empty()) {
      const int i = queue_.front();
      queue_.pop_front();
      if (cycle_[i] == pass_) continue;
      cycle_[i] = pass_;
      for (int j = 0; j < k_; ++j) {
        const int v = links_[static_cast<size_t>(i) * k_ + j].nbr;
        if (cycle_[v] != pass_) queue_.push_back(v);
      }
      if (i != start || adjust_start) {
        double err = 0.0;
        steps += adjust_point(i, err);
        total += err;
      }
    }
  }

  int n_, d_, k_;
  std::vector<double> y_;  // preserved dimensions, one point per row
  std::vector<Link> links_;
  std::vector<int> cycle_;  // iteration in which each point was last visited
  std::deque<int> queue_;
  double rate_, scale_, scale2_;
  double ave_dist_, inv_ave_, eta_;
  int pass_;
};

}  // namespace

// Fits Manifold Sculpting, `nstart` times.
//
// `x` holds the (centred) data, `y0` the starting coordinates of the preserved
// dimensions and `resid` the scaled dimensions (it may have no columns), all
// with one row per point. `nn` and `nd` are the neighbours (1-based) and their
// distances in `x`.
// [[Rcpp::export]]
Rcpp::List ms_fit_cpp(const arma::mat& x, const arma::mat& y0, const arma::mat& resid,
                      const Rcpp::IntegerMatrix& nn, const Rcpp::NumericMatrix& nd,
                      double squishing_rate, int max_iter, int patience,
                      double burn_in_scale, int nstart, bool verbose) {
  const int n = static_cast<int>(x.n_rows);
  const int d = static_cast<int>(y0.n_cols);
  const int k = nn.ncol();
  const int r = static_cast<int>(resid.n_cols);

  const std::vector<int> idx = flat_index(nn);
  const std::vector<double> dist = flat_distance(nd);
  const arma::mat xt = x.t();
  const arma::mat rt = resid.t();

  // Step 1-2: the original relationships, and the contribution of the scaled
  // dimensions to the distances and angles.
  std::vector<Link> links(static_cast<size_t>(n) * k);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < k; ++j) {
      const size_t ij = static_cast<size_t>(i) * k + j;
      Link& l = links[ij];
      l.nbr = idx[ij];
      l.dist0 = dist[ij];
      double theta = 0.0;
      l.slot = most_colinear(xt.colptr(i), l.nbr, dist[ij] * dist[ij], xt, idx, dist,
                             k, i, theta);
      l.theta0 = target_angle(theta, d);
      l.jsq = 0.0;
      l.jdot = 0.0;
      if (r > 0) {
        const double* ri = rt.colptr(i);
        const double* rv = rt.colptr(l.nbr);
        const double* rm = l.slot >= 0 ? rt.colptr(idx[static_cast<size_t>(l.nbr) * k + l.slot]) : rv;
        for (int t = 0; t < r; ++t) {
          const double a = ri[t] - rv[t];
          l.jsq += a * a;
          l.jdot += a * (rm[t] - rv[t]);
        }
      }
    }
    if (i % 256 == 0) Rcpp::checkUserInterrupt();
  }

  std::vector<double> y(static_cast<size_t>(n) * d);
  for (int i = 0; i < n; ++i) {
    for (int t = 0; t < d; ++t) y[static_cast<size_t>(i) * d + t] = y0(i, t);
  }

  // Each run starts from the same configuration, and the runs differ by their
  // random traversals. The run with the lowest final error is kept.
  Rcpp::NumericVector final_error(nstart);
  arma::mat embedding(n, d);
  std::vector<double> trace;
  int iter = 0;
  bool converged = false;
  double best_final = std::numeric_limits<double>::infinity();
  double ave_dist = 0.0;

  for (int run = 0; run < nstart; ++run) {
    Sculptor sculptor(n, d, k, y, links, squishing_rate);
    ave_dist = sculptor.ave_dist();

    // Burn-in: squish until the scaled dimensions have nearly vanished, then
    // carry on until the error has not improved for `patience` iterations.
    std::vector<double> run_trace;
    run_trace.reserve(std::min(max_iter, 10000));
    int run_iter = 0;
    bool burnt_in = false;
    int since_best = 0;
    double best = std::numeric_limits<double>::infinity();
    while (run_iter < max_iter) {
      const double err = sculptor.pass();
      run_trace.push_back(err);
      ++run_iter;
      if (!burnt_in) {
        burnt_in = sculptor.scale() <= burn_in_scale;
      } else {
        if (err < best) {
          best = err;
          since_best = 0;
        } else if (++since_best >= patience) {
          break;
        }
      }
      if (verbose && (run_iter % 50 == 0)) {
        if (nstart > 1) Rcpp::Rcout << "run " << run + 1 << ", ";
        Rcpp::Rcout << "iteration " << run_iter << ": error " << err
                    << ", scale " << sculptor.scale()
                    << ", learning rate " << sculptor.learning_rate() << "\n";
      }
      if (run_iter % 10 == 0) Rcpp::checkUserInterrupt();
    }
    sculptor.center();

    final_error[run] = run_trace.back();
    if (run_trace.back() < best_final) {
      best_final = run_trace.back();
      const std::vector<double>& yf = sculptor.y();
      for (int i = 0; i < n; ++i) {
        for (int t = 0; t < d; ++t) embedding(i, t) = yf[static_cast<size_t>(i) * d + t];
      }
      trace.swap(run_trace);
      iter = run_iter;
      converged = burnt_in && since_best >= patience;
    }
  }

  return Rcpp::List::create(
      Rcpp::Named("embedding") = embedding,
      Rcpp::Named("error") = Rcpp::wrap(trace),
      Rcpp::Named("iterations") = iter,
      Rcpp::Named("converged") = converged,
      Rcpp::Named("final_error") = final_error,
      Rcpp::Named("ave_dist") = ave_dist);
}

namespace {

// Error of a new point at `yq`, against fixed training points.
double new_point_error(const double* yq, const arma::mat& yt, const int* nbr,
                       const double* dist0, const int* mpt, const double* theta0,
                       int k, double inv_ave) {
  const int d = static_cast<int>(yt.n_rows);
  double err = 0.0;
  for (int j = 0; j < k; ++j) {
    const double* yv = yt.colptr(nbr[j]);
    double sq_a = 0.0;
    for (int t = 0; t < d; ++t) {
      const double a = yq[t] - yv[t];
      sq_a += a * a;
    }
    const double dd = (std::sqrt(sq_a) - dist0[j]) * inv_ave;
    double dth = 0.0;
    if (mpt[j] >= 0) {
      const double* ym = yt.colptr(mpt[j]);
      double dot = 0.0;
      double sq_b = 0.0;
      for (int t = 0; t < d; ++t) {
        const double a = yq[t] - yv[t];
        const double b = ym[t] - yv[t];
        dot += a * b;
        sq_b += b * b;
      }
      dth = std::min(0.0, angle01(dot, sq_a, sq_b) - theta0[j]);
    }
    err += dd * dd + dth * dth;
  }
  return err;
}

}  // namespace

// Embeds new points.
//
// Each point starts at the combination of the embeddings of its training
// neighbours given by its locally linear reconstruction weights (Roweis and
// Saul, 2000). With `refine`, it is then moved, by pattern search, to
// minimise the Manifold Sculpting error against its neighbours, whose
// embeddings are held fixed. A point that coincides with a training point
// takes its embedding.
// [[Rcpp::export]]
arma::mat ms_project_cpp(const arma::mat& x, const arma::mat& y,
                         const Rcpp::IntegerMatrix& nn, const Rcpp::NumericMatrix& nd,
                         const arma::mat& x_new, const Rcpp::IntegerMatrix& nn_new,
                         const Rcpp::NumericMatrix& nd_new, double ave_dist, bool refine) {
  const int m = static_cast<int>(x_new.n_rows);
  const int d = static_cast<int>(y.n_cols);
  const int k_train = nn.ncol();
  const int k = nn_new.ncol();

  const std::vector<int> idx = flat_index(nn);
  const std::vector<double> dist = flat_distance(nd);
  const arma::mat xt = x.t();
  const arma::mat yt = y.t();
  const arma::mat qt = x_new.t();
  const double inv_ave = 1.0 / std::max(ave_dist, 1e-300);
  const double tol = 1e-10 * ave_dist;

  arma::mat out(m, d);
  std::vector<int> nbr(k), mpt(k);
  std::vector<double> dist0(k), theta0(k);
  arma::mat z(xt.n_rows, k);
  arma::vec yq(d);

  for (int q = 0; q < m; ++q) {
    for (int j = 0; j < k; ++j) {
      nbr[j] = nn_new(q, j) - 1;
      dist0[j] = nd_new(q, j);
    }
    if (dist0[0] <= tol) {
      out.row(q) = y.row(nbr[0]);
      continue;
    }

    // Locally linear reconstruction weights, regularised as in LLE.
    const double* xq = qt.colptr(q);
    for (int j = 0; j < k; ++j) {
      const double* xv = xt.colptr(nbr[j]);
      for (arma::uword t = 0; t < xt.n_rows; ++t) z(t, j) = xv[t] - xq[t];
    }
    arma::mat gram = z.t() * z;
    const double tr = arma::trace(gram);
    gram.diag() += tr > 0.0 ? 1e-3 * tr : 1e-3;
    arma::vec w;
    bool ok = arma::solve(w, gram, arma::ones<arma::vec>(k),
                          arma::solve_opts::likely_sympd + arma::solve_opts::no_approx);
    if (!ok || !w.is_finite() || std::abs(arma::accu(w)) < 1e-12) {
      w = 1.0 / arma::vec(dist0);  // inverse distance weights
    }
    w /= arma::accu(w);
    yq.zeros();
    for (int j = 0; j < k; ++j) yq += w[j] * yt.col(nbr[j]);

    if (refine) {
      for (int j = 0; j < k; ++j) {
        double theta = 0.0;
        const int slot = most_colinear(xq, nbr[j], dist0[j] * dist0[j], xt, idx, dist,
                                       k_train, -1, theta);
        mpt[j] = slot >= 0 ? idx[static_cast<size_t>(nbr[j]) * k_train + slot] : -1;
        theta0[j] = target_angle(theta, d);
      }

      double* yp = yq.memptr();
      double base = new_point_error(yp, yt, nbr.data(), dist0.data(), mpt.data(),
                                    theta0.data(), k, inv_ave);
      double step = 0.25 * dist0[0];
      const double min_step = 1e-4 * step;
      while (step > min_step) {
        bool progress = true;
        for (int sweep = 0; progress && sweep < 100; ++sweep) {
          progress = false;
          for (int t = 0; t < d; ++t) {
            const double v = yp[t];
            yp[t] = v + step;
            double e = new_point_error(yp, yt, nbr.data(), dist0.data(), mpt.data(),
                                       theta0.data(), k, inv_ave);
            if (e >= base) {
              yp[t] = v - step;
              e = new_point_error(yp, yt, nbr.data(), dist0.data(), mpt.data(),
                                  theta0.data(), k, inv_ave);
            }
            if (e >= base) {
              yp[t] = v;
            } else {
              base = e;
              progress = true;
            }
          }
        }
        step *= 0.5;
      }
    }
    out.row(q) = yq.t();
    if (q % 256 == 0) Rcpp::checkUserInterrupt();
  }
  return out;
}
