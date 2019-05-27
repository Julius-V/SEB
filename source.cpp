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
List SEB(mat D, uvec times, const bool intervals = true, const std::string order = "successive") {
  const int nr = D.n_rows;
  const int nc = D.n_cols;
  if (times.n_elem != nc)
    stop("The number of specified cutting times has to be the same as the number of columns.");
  if (order != "successive" & order != "random") {
    stop("The order must be either successive or random.");
  } else if (order == "random") {
    const uvec o = shuffle(regspace<uvec>(0, nc - 1));
    times = times(o);
    D = D.cols(o);
  }
  const int nints = times[0];
  const int ncells = std::accumulate(times.begin(), times.end(), 1, std::multiplies<int>());
  const int ncells0 = ncells / nints;
  if (ncells > nr)
    stop("There are more cells than observations.");
  const vec D0 = D.col(0);
  const uvec idx = arma::floor((nr * regspace<uvec>(1, nints)) / nints) - 1;
  vec breaks = sort(D0);
  breaks = breaks(idx);
  breaks(nints - 1) = R_PosInf;
  const uvec labelsVec = fastCut(D0, breaks);
  if (nc > 1) {
    uvec labelsMat(nr);
    const mat D1 = D.cols(1, nc - 1);
    mat Ints(ncells, nc * 2);
    for (int i = 1; i <= nints; i++) {
      const mat subD = D1.rows(find(labelsVec == i));
      const List subSEB = SEB(subD, times.subvec(1, nc - 1), intervals);
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
      mat Ints0(nints, 2);
      Ints0(0, 0) = R_NegInf;
      if(nints > 1)
        Ints0.submat(1, 0, nints - 1, 0) = breaks.subvec(0, nints - 2);
      Ints0.col(1) = breaks;
      return List::create(_["labels"] = as<std::vector<int> >(wrap(labelsVec)), _["intervals"] = Ints0);
    } else {
      return List::create(_["labels"] = as<std::vector<int> >(wrap(labelsVec)));  
    }
  }
}