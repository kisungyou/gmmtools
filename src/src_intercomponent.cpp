/* INTER-COMPONENT DISTANCE COMPUTATION
 * (01) interdist_bhat : Bhattacharyya Distance
 */

#include "RcppArmadillo.h"
#include "utilities.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// (01) Bhattacharyya Distance -------------------------------------------------
// [[Rcpp::export]]
arma::mat interdist_bhat(arma::mat &par_means, arma::cube &par_covs){
  // parameters
  int n = par_means.n_rows;
  int p = par_means.n_cols;
  
  // pre-compute
  arma::vec logdets(n,fill::zeros);
  for (int i=0; i<n; i++){
    logdets(i) = arma::log_det_sympd(par_covs.slice(i));
  }
  
  // main computation
  arma::colvec mudiff(p,fill::zeros);
  arma::mat meancov(p,p,fill::zeros);
  
  double term1 = 0.0;
  double term2 = 0.0;
  
  arma::mat output(n,n,fill::zeros);
  for (int i=0; i<(n-1); i++){
    for (int j=(i+1); j<n; j++){
      // setup
      mudiff  = arma::trans(par_means.row(i)-par_means.row(j));
      meancov = (par_covs.slice(i)+par_covs.slice(j))/2.0;
    
      // compute
      term1 = arma::dot(arma::solve(meancov, mudiff), mudiff)/8.0;
      term2 = 0.5*arma::log_det_sympd(meancov);
      
      // fill-in
      output(i,j) = term1 + term2 - (0.25*(logdets(i) + logdets(j)));
      output(j,i) = output(i,j);
    }
  }
  return(output);
}