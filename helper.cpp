#include <cmath>
#include <Rmath.h>
#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
mat mMmu_x(mat beta_x, 
           IntegerVector idx,
           NumericVector mMmu,
           bool add){
  IntegerVector mDims = mMmu.attr("dim");
  mat a_mMmu(mMmu.begin(), mDims[0], mDims[1]);
  // update mMmu
  for (int j = 0; j < mDims[1]; j++){
    if (add){
      a_mMmu.col(j) += beta_x.col( idx[j]-1 );
    }else{
      a_mMmu.col(j) -= beta_x.col( idx[j]-1 );
    }
  }
  return(a_mMmu);
}

// [[Rcpp::export]]
mat mMmu_M(NumericVector beta_M,
           NumericVector mMmu,
           bool add){
  IntegerVector mDims = mMmu.attr("dim");
  mat a_mMmu(mMmu.begin(), mDims[0], mDims[1]);
  vec a_beta_M(beta_M.begin(), mDims[0], false);
  // update mMmu
  for (int i = 0; i < mDims[1]; i++){
    if (add){
      a_mMmu.col(i) += a_beta_M;
    }else{
      a_mMmu.col(i) -= a_beta_M;
    }
  }
  return(a_mMmu);
}