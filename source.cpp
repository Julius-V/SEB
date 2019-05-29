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
List SEB(mat D, uvec nints, const bool intervals = true, const std::string order = "successive") {
  const int nr = D.n_rows;
  const int nc = D.n_cols;
  if (nints.n_elem != nc)
    stop("'nints' has to have as many elements as 'D' has has columns.");
  if (order != "successive" & order != "random") {
    stop("'order' must be either successive or random.");
  } else if (order == "random") {
    const uvec o = shuffle(regspace<uvec>(0, nc - 1));
    nints = nints(o);
    D = D.cols(o);
  }
  const int nints0 = nints[0];
  const int ncells = std::accumulate(nints.begin(), nints.end(), 1, std::multiplies<int>());
  const int ncells0 = ncells / nints0;
  if (ncells > nr)
    stop("The total number of cells cannot exceed the number of observations.");
  const vec D0 = D.col(0);
  const uvec idx = arma::floor((nr * regspace<uvec>(1, nints0)) / nints0) - 1;
  vec breaks = sort(D0);
  breaks = breaks(idx);
  breaks(nints0 - 1) = R_PosInf;
  const uvec labelsVec = fastCut(D0, breaks);
  if (nc > 1) {
    uvec labelsMat(nr);
    const mat D1 = D.cols(1, nc - 1);
    mat Ints(ncells, nc * 2);
    for (int i = 1; i <= nints0; i++) {
      const mat subD = D1.rows(find(labelsVec == i));
      const List subSEB = SEB(subD, nints.subvec(1, nc - 1), intervals);
      const uvec subLabels = subSEB[0];
      int k = 0;
      for (int j = 0; j < nr; j++) {
        if (labelsVec[j] == i) {
          labelsMat[j] = subLabels[k] + ncells0 * (i - 1);
          k++;
        }
      }
      if (intervals) {
        const mat subInts = subSEB["intervals"];
        Ints.submat(ncells0 * (i - 1), 0, ncells0 * i - 1, 0).fill((i == 1) ? R_NegInf : breaks[i - 2]);
        Ints.submat(ncells0 * (i - 1), 1, ncells0 * i - 1, 1).fill(breaks[i - 1]);
        Ints.submat(ncells0 * (i - 1), 2, ncells0 * i - 1, 2 * nc - 1) = subInts;
      }
    }
    if (intervals) {
      return List::create(_["labels"] = as<std::vector<int> >(wrap(labelsMat)), _["intervals"] = Ints);
    } else {
      return List::create(_["labels"] = as<std::vector<int> >(wrap(labelsMat)));  
    }
  } else {
    if (intervals) {
      mat Ints0(nints0, 2);
      Ints0(0, 0) = R_NegInf;
      if(nints0 > 1)
        Ints0.submat(1, 0, nints0 - 1, 0) = breaks.subvec(0, nints0 - 2);
      Ints0.col(1) = breaks;
      return List::create(_["labels"] = as<std::vector<int> >(wrap(labelsVec)), _["intervals"] = Ints0);
    } else {
      return List::create(_["labels"] = as<std::vector<int> >(wrap(labelsVec)));  
    }
  }
}