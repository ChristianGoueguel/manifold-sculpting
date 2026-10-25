// Exact k-nearest neighbour search (Euclidean distance).
//
// Two strategies are available: a k-d tree, which is fast when the data have
// few columns, and a blocked brute-force search, whose cross-products are
// computed by BLAS and which is the better choice for many columns.

#include <RcppArmadillo.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <utility>
#include <vector>

namespace {

// A candidate neighbour: (squared distance, index). Pairs are compared
// lexicographically, so ties in distance are broken by the smaller index and
// both search strategies return the same neighbours.
typedef std::pair<double, int> Candidate;

inline double sq_dist(const double* a, const double* b, int p) {
  double s = 0.0;
  for (int t = 0; t < p; ++t) {
    const double d = a[t] - b[t];
    s += d * d;
  }
  return s;
}

// Bounded max-heap that keeps the k best candidates seen so far.
class KBest {
 public:
  explicit KBest(int k) : k_(k) { heap_.reserve(k); }

  void clear() { heap_.clear(); }

  double worst() const {
    if (static_cast<int>(heap_.size()) < k_) {
      return std::numeric_limits<double>::infinity();
    }
    return heap_.front().first;
  }

  void push(double d2, int j) {
    const Candidate c(d2, j);
    if (static_cast<int>(heap_.size()) < k_) {
      heap_.push_back(c);
      std::push_heap(heap_.begin(), heap_.end());
    } else if (c < heap_.front()) {
      std::pop_heap(heap_.begin(), heap_.end());
      heap_.back() = c;
      std::push_heap(heap_.begin(), heap_.end());
    }
  }

  // The candidates, closest first.
  const std::vector<Candidate>& sorted() {
    std::sort_heap(heap_.begin(), heap_.end());
    return heap_;
  }

 private:
  int k_;
  std::vector<Candidate> heap_;
};

// k-d tree over the columns of a p x n matrix (one point per column), split
// at the median of the coordinate with the largest spread.
class KdTree {
 public:
  KdTree(const double* pts, int n, int p, int leaf_size)
      : pts_(pts), p_(p), leaf_size_(leaf_size), index_(n) {
    std::iota(index_.begin(), index_.end(), 0);
    nodes_.reserve(2 * (n / leaf_size + 1));
    build(0, n);
  }

  // Adds the neighbours of q to `best`, skipping the point `exclude`.
  void query(const double* q, int exclude, KBest& best) const {
    search(0, q, exclude, best);
  }

 private:
  struct Node {
    int lo, hi;    // range of index_ covered by the node
    int dim;       // split coordinate, -1 for a leaf
    double split;  // split value
    int left, right;
  };

  double coord(int i, int d) const { return pts_[static_cast<size_t>(i) * p_ + d]; }

  int build(int lo, int hi) {
    const int id = static_cast<int>(nodes_.size());
    nodes_.push_back(Node{lo, hi, -1, 0.0, -1, -1});
    if (hi - lo <= leaf_size_) {
      return id;
    }

    int dim = -1;
    double spread = 0.0;
    for (int d = 0; d < p_; ++d) {
      double lo_v = coord(index_[lo], d);
      double hi_v = lo_v;
      for (int t = lo + 1; t < hi; ++t) {
        const double v = coord(index_[t], d);
        lo_v = std::min(lo_v, v);
        hi_v = std::max(hi_v, v);
      }
      if (hi_v - lo_v > spread) {
        spread = hi_v - lo_v;
        dim = d;
      }
    }
    if (dim < 0) {
      return id;  // all the points are identical
    }

    const int mid = lo + (hi - lo) / 2;
    std::nth_element(index_.begin() + lo, index_.begin() + mid, index_.begin() + hi,
                     [this, dim](int a, int b) { return coord(a, dim) < coord(b, dim); });
    const double split = coord(index_[mid], dim);
    const int left = build(lo, mid);
    const int right = build(mid, hi);

    Node& node = nodes_[id];
    node.dim = dim;
    node.split = split;
    node.left = left;
    node.right = right;
    return id;
  }

  void search(int id, const double* q, int exclude, KBest& best) const {
    const Node& node = nodes_[id];
    if (node.dim < 0) {
      for (int t = node.lo; t < node.hi; ++t) {
        const int j = index_[t];
        if (j != exclude) {
          best.push(sq_dist(q, pts_ + static_cast<size_t>(j) * p_, p_), j);
        }
      }
      return;
    }
    // Points left of the split have coordinates <= split, points right of it
    // have coordinates >= split.
    const double diff = q[node.dim] - node.split;
    const int near = diff <= 0.0 ? node.left : node.right;
    const int far = diff <= 0.0 ? node.right : node.left;
    search(near, q, exclude, best);
    if (diff * diff <= best.worst()) {
      search(far, q, exclude, best);
    }
  }

  const double* pts_;
  int p_;
  int leaf_size_;
  std::vector<int> index_;
  std::vector<Node> nodes_;
};

}  // namespace

// k nearest neighbours in `ref` of each row of `query`.
//
// Returns the 1-based indices of the neighbours and their Euclidean
// distances, as m x k matrices sorted by distance. With `exclude_self`,
// `query` must be `ref` and each point is left out of its own neighbours.
// [[Rcpp::export]]
Rcpp::List knn_cpp(const arma::mat& ref, const arma::mat& query, int k,
                   bool exclude_self, bool use_tree) {
  const int n = static_cast<int>(ref.n_rows);
  const int m = static_cast<int>(query.n_rows);
  const int p = static_cast<int>(ref.n_cols);
  if (static_cast<int>(query.n_cols) != p) {
    Rcpp::stop("`ref` and `query` must have the same number of columns.");
  }
  if (k < 1 || k > n - (exclude_self ? 1 : 0)) {
    Rcpp::stop("`k` must be between 1 and the number of reference points.");
  }

  Rcpp::IntegerMatrix index(m, k);
  Rcpp::NumericMatrix distance(m, k);

  // One point per column, so that the coordinates of a point are contiguous.
  const arma::mat ref_t = ref.t();
  const arma::mat query_t = query.t();

  if (use_tree) {
    KdTree tree(ref_t.memptr(), n, p, 16);
    KBest best(k);
    for (int i = 0; i < m; ++i) {
      best.clear();
      tree.query(query_t.colptr(i), exclude_self ? i : -1, best);
      const std::vector<Candidate>& nb = best.sorted();
      for (int j = 0; j < k; ++j) {
        index(i, j) = nb[j].second + 1;
        distance(i, j) = std::sqrt(nb[j].first);
      }
      if (i % 1024 == 0) Rcpp::checkUserInterrupt();
    }
  } else {
    // Squared distances, up to the squared norm of the query (a constant for
    // each query), are |r|^2 - 2 r.q. A few extra candidates are kept and
    // ranked again by their exact distances, which corrects the rounding
    // errors of this expansion.
    const arma::vec ref_norm = arma::sum(arma::square(ref), 1);
    const int extra = std::min(k + 8, n - (exclude_self ? 1 : 0));
    const int block = std::max(1, std::min(m, static_cast<int>(8388608 / std::max(n, 1))));
    KBest approx(extra);
    KBest best(k);

    for (int b0 = 0; b0 < m; b0 += block) {
      const int b1 = std::min(m, b0 + block);
      const arma::mat cross = ref * query.rows(b0, b1 - 1).t();  // n x (b1 - b0)
      for (int c = 0; c < b1 - b0; ++c) {
        const int i = b0 + c;
        const double* g = cross.colptr(c);
        approx.clear();
        for (int j = 0; j < n; ++j) {
          const double score = ref_norm[j] - 2.0 * g[j];
          if (score <= approx.worst() && !(exclude_self && j == i)) {
            approx.push(score, j);
          }
        }
        best.clear();
        const double* q = query_t.colptr(i);
        for (const Candidate& cand : approx.sorted()) {
          best.push(sq_dist(q, ref_t.colptr(cand.second), p), cand.second);
        }
        const std::vector<Candidate>& nb = best.sorted();
        for (int j = 0; j < k; ++j) {
          index(i, j) = nb[j].second + 1;
          distance(i, j) = std::sqrt(nb[j].first);
        }
      }
      Rcpp::checkUserInterrupt();
    }
  }

  return Rcpp::List::create(Rcpp::Named("index") = index,
                            Rcpp::Named("distance") = distance);
}
