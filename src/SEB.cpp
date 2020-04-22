#define RCPP_ARMADILLO_RETURN_COLVEC_AS_VECTOR
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

uvec fastCut(vec x, vec breaks) {
  const int n = x.n_elem;
  uvec labels(n);
  for (int i = 0; i < n; i++) {
    int j = 0;
    while (x[i] > breaks[j]) {
      j++;
    }
    labels[i] = j + 1;
  }
  return labels;
}

// [[Rcpp::export]]
List SEB(arma::mat X, arma::uvec nints, const bool intervals = true, const std::string order = "successive") {
  const int nr = X.n_rows;
  const int nc = X.n_cols;
  if (order == "random") {
    const uvec o = shuffle(regspace<uvec>(0, nc - 1));
    nints = nints(o);
    X = X.cols(o);
  }
  const int nints0 = nints[0];
  const int nblocks = std::accumulate(nints.begin(), nints.end(), 1, std::multiplies<int>());
  const int nblocks1 = nblocks / nints0;
  const vec X0 = X.col(0);
  const uvec idx = arma::floor((nr * regspace<uvec>(1, nints0)) / nints0) - 1;
  vec breaks = sort(X0);
  breaks = breaks(idx);
  breaks(nints0 - 1) = R_PosInf;
  const uvec labelsVec = fastCut(X0, breaks);
  if (nc > 1) {
    uvec labelsMat(nr);
    const mat X1 = X.cols(1, nc - 1);
    mat Ints(nblocks, nc * 2);
    for (int i = 1; i <= nints0; i++) {
      const mat subX = X1.rows(find(labelsVec == i));
      const List subSEB = SEB(subX, nints.subvec(1, nc - 1), intervals);
      const uvec subLabels = subSEB[0];
      int k = 0;
      for (int j = 0; j < nr; j++) {
        if (labelsVec[j] == i) {
          labelsMat[j] = subLabels[k] + nblocks1 * (i - 1);
          k++;
        }
      }
      if (intervals) {
        const mat subInts = subSEB["intervals"];
        Ints.submat(nblocks1 * (i - 1), 0, nblocks1 * i - 1, 0).fill((i == 1) ? R_NegInf : breaks[i - 2]);
        Ints.submat(nblocks1 * (i - 1), 1, nblocks1 * i - 1, 1).fill(breaks[i - 1]);
        Ints.submat(nblocks1 * (i - 1), 2, nblocks1 * i - 1, 2 * nc - 1) = subInts;
      }
    }
    if (intervals) {
      return List::create(_["labels"] = as<IntegerVector>(wrap(labelsMat)), _["intervals"] = Ints);
    } else {
      return List::create(_["labels"] = as<IntegerVector>(wrap(labelsMat)));  
    }
  } else {
    if (intervals) {
      mat Ints0(nints0, 2);
      Ints0(0, 0) = R_NegInf;
      if(nints0 > 1)
        Ints0.submat(1, 0, nints0 - 1, 0) = breaks.subvec(0, nints0 - 2);
      Ints0.col(1) = breaks;
      return List::create(_["labels"] = as<IntegerVector>(wrap(labelsVec)), _["intervals"] = Ints0);
    } else {
      return List::create(_["labels"] = as<IntegerVector>(wrap(labelsVec)));  
    }
  }
}