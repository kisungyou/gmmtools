#include "RcppArmadillo.h"
#include "utilities.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

/* UTILITIES FOR GMM
 * (01) gmm_density     : evaluate density of GMM
 * (02) gmm_sample      : sample from GMM
 * (03) gmm_compdist_W2 : between-component distance using 2-Wasserstein distance
 */

// (01) gmm_density ============================================================
// [[Rcpp::export]]
arma::vec gmm_density(arma::mat &coords, arma::vec &weight, arma::mat &mean, arma::cube &variance){
  // parameters
  int N = coords.n_rows;
  int K = weight.n_elem;
  
  // evaluate per Gaussian + weight
  arma::vec myweight = weight/arma::accu(weight);
  arma::mat eval_class(N,K,fill::zeros);
  for (int k=0; k<K; k++){
    eval_class.col(k) = eval_gaussian_multiple(coords, mean.row(k), variance.slice(k), false)*myweight(k);
  }
  
  // finalize
  arma::vec output(N,fill::zeros);
  for (int n=0; n<N; n++){
    output(n) = arma::accu(eval_class.row(n));
  }
  return(output);
}
// (02) gmm_sample  : sample from GMM ==========================================
// [[Rcpp::export]]
arma::mat gmm_sample(int n, arma::colvec oldweight, arma::mat oldmeans, arma::cube oldcovs){
  // model sizes
  int k = oldcovs.n_slices;
  int p = oldcovs.n_cols;
  
  arma::gmm_full model;
  model.reset(p, k);
  model.set_means(arma::trans(oldmeans)); // column centroids
  model.set_fcovs(oldcovs);
  model.set_hefts(arma::trans(oldweight));
  
  arma::mat output = arma::trans(model.generate(n));
  return(output);
}


// (03) gmm_compdist_W2 ========================================================
// [[Rcpp::export]]
arma::mat gmm_compdist_W2(arma::mat& means, arma::cube& vars){
  // parameters
  int N = means.n_rows; // number of observations
  int p = means.n_cols; // dimension
  
  // preliminary computation
  arma::cube varsqrt(p,p,N,fill::zeros);
  for (int n=0; n<N; n++){
    varsqrt.slice(n) = arma::sqrtmat_sympd(vars.slice(n));
  }
  
  // main distance computation
  arma::mat output(N,N,fill::zeros);
  for (int i=0; i<(N-1); i++){
    for (int j=(i+1); j<N; j++){
      output(i,j) = gauss2dist_wass2(means.row(i), vars.slice(i), means.row(j), vars.slice(j), varsqrt.slice(j));
      output(j,i) = output(i,j);
    }
  }
  return(output);
}
