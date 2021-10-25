#include "RcppArmadillo.h"
#include "utilities.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

/* DISSIMILARITY MEASURES
 * (01) cpp_gmmdist_hd   : Kernel/Hilbert Embedding
 * (02) cpp_gmmdist_l2   : Standard L2 norm
 * 
 * 
 */


// (01) cpp_gmmdist_hd   : Kernel/Hilbert Embedding ============================
// evaluate the inner product of the kernel \int k(x,y)dP(x)dP(y)
double cpp_gmmdist_hd_inner(arma::rowvec m1, arma::mat s1, arma::rowvec m2, arma::mat s2, double theta){
  // int d = m1.n_elem; 
  // double dd = static_cast<double>(d);
  // 
  // arma::mat sig_theta  = s1+s2+((theta*theta)*(arma::eye<arma::mat>(d,d)));
  // double tmpval   = eval_gaussian_single(m1, m2, sig_theta, true);
  // double adjuster = (dd/2.0)*std::log(2.0*arma::datum::pi) + dd*std::log(theta);
  // double output   = std::exp(tmpval + adjuster);
  // return(output);
  int d = m1.n_elem; double dd = static_cast<double>(d);
  double adjpi  = std::pow((2.0*arma::datum::pi), dd/2.0);
  double thetad = std::pow(theta, dd);
  
  arma::mat tmpcov = s1 + s2 + (theta*theta)*arma::eye<arma::mat>(d,d);
  double output = eval_gaussian_single(m1,m2,tmpcov,false)*adjpi*thetad;
  return(output);
}

// [[Rcpp::export]]
double cpp_gmmdist_hd(arma::vec weight1, arma::mat mean1, arma::cube variance1,
                      arma::vec weight2, arma::mat mean2, arma::cube variance2, double theta){
  // parameters
  int N = weight1.n_elem;
  int M = weight2.n_elem;
  
  // compute three terms
  double term1 = 0.0;
  double term2 = 0.0;
  double term3 = 0.0;
  
  for (int i=0; i<(N-1); i++){
    for (int j=(i+1); j<N; j++){
      term1 += 2.0*cpp_gmmdist_hd_inner(mean1.row(i), variance1.slice(i), mean1.row(j), variance1.slice(j), theta)*weight1(i)*weight1(j);
    }
  }
  for (int i=0; i<N; i++){
    term1 += cpp_gmmdist_hd_inner(mean1.row(i), variance1.slice(i), mean1.row(i), variance1.slice(i), theta)*(weight1(i)*weight1(i));
  }
  for (int i=0; i<(M-1); i++){
    for (int j=(i+1); j<M; j++){
      term2 += 2.0*cpp_gmmdist_hd_inner(mean2.row(i), variance2.slice(i), mean2.row(j), variance2.slice(j), theta)*weight2(i)*weight2(j);
    }
  }
  for (int i=0; i<M; i++){
    term2 += cpp_gmmdist_hd_inner(mean2.row(i), variance2.slice(i), mean2.row(i), variance2.slice(i),theta)*(weight2(i)*weight2(i));
  }
  for (int i=0; i<N; i++){
    for (int j=0; j<M; j++){
      term3 += cpp_gmmdist_hd_inner(mean1.row(i), variance1.slice(i), mean2.row(j), variance2.slice(j), theta)*weight1(i)*weight2(j);
    }
  }
  double output = std::sqrt(term1+term2-(2.0*term3));
  return(output);
}

// (02) cpp_gmmdist_l2   : Standard L2 norm ====================================
// [[Rcpp::export]]
double cpp_gmmdist_l2(arma::vec weight1, arma::mat mean1, arma::cube variance1,
                      arma::vec weight2, arma::mat mean2, arma::cube variance2){
  
  // PARAMETERS
  int N1 = weight1.n_elem;
  int N2 = weight2.n_elem;
  // int p  = variance1.n_rows;

  // COMPUTE THREE TERMS
  // GMM No.1
  double term1 = 0.0;
  for (int i=0; i<(N1-1); i++){
    for (int j=(i+1); j<N1; j++){
      term1 += 2.0*eval_gaussian_single(mean1.row(i), mean1.row(j), variance1.slice(i)+variance1.slice(j), false)*weight1(i)*weight1(j);
    }
  }
  for (int i=0; i<N1; i++){
    term1 += eval_gaussian_single(mean1.row(i), mean1.row(i), 2.0*variance1.slice(i), false)*(weight1(i)*weight1(i));
  }
  // GMM No.2
  double term2 = 0.0;
  for (int i=0; i<(N2-1); i++){
    for (int j=(i+1); j<N2; j++){
      term2 += 2.0*eval_gaussian_single(mean2.row(i), mean2.row(j), variance2.slice(i)+variance2.slice(j), false)*weight2(i)*weight2(j);
    }
  }
  for (int i=0; i<N2; i++){
    term2 += eval_gaussian_single(mean2.row(i), mean2.row(i), 2.0*variance2.slice(i), false)*(weight2(i)*weight2(i));
  }
  // CROSS TERMS
  double term3 = 0.0;
  for (int i=0; i<N1; i++){
    for (int j=0; j<N2; j++){
      term3 += eval_gaussian_single(mean1.row(i), mean2.row(j), variance1.slice(i)+variance2.slice(j), false)*weight1(i)*weight2(j);
    }
  }
  
  double output = std::sqrt(term1+term2-(2.0*term3));
  return(output);
}