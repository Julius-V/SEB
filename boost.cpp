#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
NumericVector stlSort(NumericVector x) {
  NumericVector y = clone(x);
  std::sort(y.begin(), y.end());
  return y;
}

// [[Rcpp::export]]
NumericMatrix empMat(NumericVector V, NumericVector nus) {
  const int nr = V.size();
  const int nc = nus.size();
  NumericMatrix mat(nr, nc);
  for (int i = 0; i < nr; i++) {
    for (int j = 0; j < nc; j++) {
      mat(i, j) = (V[i] <= nus[j]) - nus[j];
    }
  }
  return mat;
}

// [[Rcpp::export]]
NumericVector armaColsums(const arma::mat& x) {
  arma::rowvec y = arma::sum(x, 0);
  return NumericVector(y.begin(), y.end());
}

// [[Rcpp::export]]
NumericVector armaOrder(arma::vec x) {
  NumericVector y = as<NumericVector>(wrap(arma::sort_index(x) + 1));
  return NumericVector(y.begin(), y.end());
}

// [[Rcpp::export]]
LogicalMatrix indMat1(NumericVector X) {
  const int n = X.size();
  LogicalMatrix mat(n, n);
  for (int j = 0; j < n; j++) {
    mat(_, j) = X <= X(j);
  }
  return mat;
}

// [[Rcpp::export]]
LogicalMatrix indMat2(NumericVector X, NumericVector C, IntegerVector idx) {
  const int nr = X.size();
  const int nc = idx.size();
  LogicalMatrix mat(nr, nc);
  LogicalVector fv (nr, false);
  for (int j = 0; j < nc; j++) {
    mat(_, j) = (idx(j) >= 0) ? (X <= C(idx(j))) : fv;
  }
  return mat;
}

// [[Rcpp::export]]
NumericMatrix hadamard(NumericMatrix A, NumericMatrix B) {
  const int nr = B.nrow();
  const int nc = B.ncol();
  NumericMatrix mat(nr, nc);
  for (int i = 0; i < nr; i++) {
    for (int j = 0; j < nc; j++) {
      mat(i, j) = A(i, j) * B(i, j);
    }
  }
  return mat;
}

// [[Rcpp::export]]
float aCKtildeU(NumericVector V, IntegerVector o) {
  const int n = V.size();
  NumericVector out(n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j <= i; j++) {
      out[i] += (V[o[j]] <= V[i]) - V[i];
    }
  }
  return max(abs(out));
}

// [[Rcpp::export]]
float aCKtilde(NumericVector V, NumericMatrix X) {
  const int n = V.size();
  const int nc = X.ncol();
  NumericVector out(n);
  int go;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      go = 0;
      for (int k = 0; k < nc; k++) {
        go +=  X(j, k) <= X(i, k);
      }
      out[i] += (go == nc) ? ((V[j] <= V[i]) - V[i]) : 0;
    }
  }
  return max(abs(out));
}

IntegerVector fastCut(NumericVector x, NumericVector breaks) {
  int n = x.size();
  IntegerVector labels(n);
  for (int i = 0; i < n; i++) {
    int j = 0;
    while (x[i] > breaks[j]) {
      j++;
    }
    labels[i] = j + 1;
  }
  return labels;
}

NumericMatrix fastSubset(NumericMatrix X, NumericVector cond) {
  mat aux(X.begin(), X.nrow(), X.ncol(), false);
  colvec condIdx(cond.begin(), cond.size(), false); 
  mat subX = aux.rows(find(condIdx == 1));
  return wrap(subX);
}

NumericMatrix fastPermute(NumericMatrix X, IntegerVector o) {
  mat matX = as<mat>(X);
  matX = matX.cols(as<uvec>(o));
  return wrap(matX);
}

// [[Rcpp::export]]
IntegerVector fastSEB(NumericMatrix D, int C, Rcpp::Nullable<int> L = R_NilValue) {
  const int nr = D.nrow();
  const int nc = D.ncol();
  NumericVector D0 = D(_, 0);
  int nints = Rf_isNull(L) ? C : Rcpp::as<int>(L);
  NumericVector breaks = stlSort(D0)[floor((nr * seq_len(nints)) / nints) - 1];
  IntegerVector labels0 = fastCut(D0, breaks);
  if (nc > 1) {
    IntegerVector labels(nr);
    NumericMatrix D1 = D(_, Range(1, nc - 1));
    for (int i = 1; i <= nints; i++) {
      NumericMatrix subD = fastSubset(D1, wrap(labels0 == i));
      IntegerVector subLabels = fastSEB(subD, C);
      int k = 0;
      for (int j = 0; j < nr; j++) {
        if (labels0[j] == i) {
          labels(j) = subLabels[k] + pow(C, nc - 1) * (i - 1);
          k++;
        }
      }
      
    }
    return labels;
  } else {
    return labels0;
  }
}

// [[Rcpp::export]]
List fastSEB2(NumericMatrix D, int C, Rcpp::Nullable<int> L = R_NilValue) {
  const int nr = D.nrow();
  const int nc = D.ncol();
  NumericVector D0 = D(_, 0);
  int nints = Rf_isNull(L) ? C : Rcpp::as<int>(L);
  if (nints * pow(C, nc - 1) > nr) {
    stop("More cells than observations.");
  }
  NumericVector breaks = stlSort(D0)[floor((nr * seq_len(nints)) / nints) - 1];
  breaks[nints - 1] = R_PosInf;
  IntegerVector labels0 = fastCut(D0, breaks);
  if (nc > 1) {
    IntegerVector labels(nr);
    NumericMatrix D1 = D(_, Range(1, nc - 1));
    mat Ints(nints * pow(C, nc - 1), nc * 2);
    for (int i = 1; i <= nints; i++) {
      NumericMatrix subD = fastSubset(D1, wrap(labels0 == i));
      List subLabels = fastSEB2(subD, C);
      IntegerVector subLabels0 = subLabels[0];
      mat subInts = subLabels["Ints"];
      int k = 0;
      for (int j = 0; j < nr; j++) {
        if (labels0[j] == i) {
          labels(j) = subLabels0[k] + pow(C, nc - 1) * (i - 1);
          k++;
        }
      }
      Ints.submat(pow(C, nc - 1) * (i - 1), 0, pow(C, nc - 1) * i - 1, 0).fill((i == 1) ? R_NegInf : breaks[i - 2]);
      Ints.submat(pow(C, nc - 1) * (i - 1), 1, pow(C, nc - 1) * i - 1, 1).fill(breaks[i - 1]);
      Ints.submat(pow(C, nc - 1) * (i - 1), 2, pow(C, nc - 1) * i - 1, 2 * nc - 1) = subInts;
      // return List::create(Ints);
    }
    if (Rf_isNull(L)) {
      return List::create(Named("labels") = labels, Named("Ints") = Ints);
    } else {
      return List::create(Named("labels0") = labels0, Named("labels") = labels, Named("Ints") = Ints);
    }
  } else {
    mat Ints0(nints, 2);
    Ints0(0, 0) = R_NegInf;
    Ints0.submat(1, 0, nints - 1, 0) = (as<vec>(breaks)).subvec(0, nints - 2);
    Ints0.col(1) = as<vec>(breaks);
    return List::create(Named("labels0") = labels0, Named("Ints") = Ints0);
  }
}

// [[Rcpp::export]]
List fastSEB3(NumericMatrix D, IntegerVector times, std::string order = "successive") {
  const int nr = D.nrow();
  const int nc = D.ncol();
  if (times.size() != nc) {
    stop("The number of specified cutting times has to be the same as the number of columns.");
  }
  if (order != "successive" & order != "random") {
    stop("The order must be either successive or random.");
  } else if (order == "random") {
    IntegerVector o = sample(nc, nc) - 1;
    times = times[o];
    D = fastPermute(D, o);
  }
  NumericVector D0 = D(_, 0);
  const int nints = times[0];
  const int ncells = std::accumulate(times.begin(), times.end(), 1, std::multiplies<int>());
  const int ncells0 = ncells / nints;
  if (ncells > nr) {
    stop("More cells than observations.");
  }
  NumericVector breaks = stlSort(D0)[Rcpp::floor((nr * seq_len(nints)) / nints) - 1];
  breaks[nints - 1] = R_PosInf;
  IntegerVector labels0 = fastCut(D0, breaks);
  if (nc > 1) {
    IntegerVector labels(nr);
    NumericMatrix D1 = D(_, Range(1, nc - 1));
    mat Ints(ncells, nc * 2);
    for (int i = 1; i <= nints; i++) {
      NumericMatrix subD = fastSubset(D1, wrap(labels0 == i));
      List subLabels = fastSEB3(subD, times[Range(1, nc - 1)]);
      IntegerVector subLabels0 = subLabels[0];
      mat subInts = subLabels["Ints"];
      int k = 0;
      for (int j = 0; j < nr; j++) {
        if (labels0[j] == i) {
          labels(j) = subLabels0[k] + ncells0 * (i - 1);
          k++;
        }
      }
      Ints.submat(ncells0 * (i - 1), 0, ncells0 * i - 1, 0).fill((i == 1) ? R_NegInf : breaks[i - 2]);
      Ints.submat(ncells0 * (i - 1), 1, ncells0 * i - 1, 1).fill(breaks[i - 1]);
      Ints.submat(ncells0 * (i - 1), 2, ncells0 * i - 1, 2 * nc - 1) = subInts;
    }
    return List::create(_["labels"] = labels, _["Ints"] = Ints);
  } else {
    mat Ints0(nints, 2);
    Ints0(0, 0) = R_NegInf;
    Ints0.submat(1, 0, nints - 1, 0) = (as<vec>(breaks)).subvec(0, nints - 2);
    Ints0.col(1) = as<vec>(breaks);
    return List::create(_["labels"] = labels0, _["Ints"] = Ints0);
  }
}