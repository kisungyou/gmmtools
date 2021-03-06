#include "RcppArmadillo.h"
#include "utilities.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

/* FITTING GMM MODELS
 * (01) gmm_armadillo : use Armadillo's native routine for Gaussian Mixture
 * (02) gmm_11R       : regularized version of Ruan et al. (2011)
 * (03) gmm_16Gfix    : use fixed weights
 * 
 * 
 * GENERIC ROUTINES
 * gmm_loglkd         : log-likelihood of the data given the model
 * gmm_standard_gamma : update Gamma    (E-STEP)
 * gmm_standard_pi    : update Pi       (M-STEP)
 * gmm_standard_mean  : update Mean     (M-STEP)
 * gmm_standard_cov   : update Variance (M-STEP)
 * gmm_skeleton       : exemplary skeletal code for GMM
 * gmm_combine_wsum   : weighted sum
 * gmm_density        : density of a gmm model
 * gmm_pdist_wass2    : compute pairwise distance between GMM components (Wass2)
 * gmm_w2barycenter   : W2 barycenter of Gaussian distributions
 */


double gmm_loglkd(arma::mat X, arma::colvec oldweight, arma::mat oldmeans, arma::cube oldcovs){
  // model sizes
  int k = oldcovs.n_slices;
  int p = oldcovs.n_cols;
  
  arma::gmm_full model;
  model.reset(p, k);
  model.set_means(arma::trans(oldmeans)); // column centroids
  model.set_fcovs(oldcovs);
  model.set_hefts(arma::trans(oldweight));
  
  double output = model.sum_log_p(arma::trans(X));
  return(output);
}
arma::mat gmm_standard_gamma(arma::mat X, arma::mat kMu, arma::cube kSig, arma::vec kPi){
  // parameters
  int N = X.n_rows;
  int K = kPi.n_elem;
  int P = X.n_cols;
  
  arma::mat probmat(N,K,fill::zeros);
  arma::rowvec parmu(P,fill::zeros);
  arma::mat    parsig(P,P,fill::zeros);
  for (int k=0;k<K;k++){
    parmu  = kMu.row(k);
    parsig = kSig.slice(k);
    probmat.col(k) = eval_gaussian_multiple(X, parmu, parsig, false)*kPi(k);
  }
  
  arma::rowvec vecK(K,fill::zeros);
  for (int n=0;n<N;n++){
    vecK = probmat.row(n);
    probmat.row(n) = vecK/arma::accu(vecK);
  }
  return(probmat);
}
arma::vec gmm_standard_pi(arma::mat Gamma){
  // parameters
  int N = Gamma.n_rows; double NN = static_cast<double>(N);
  int K = Gamma.n_cols;
  arma::vec outPI(K,fill::zeros);
  for (int k=0;k<K;k++){
    outPI(k) = arma::accu(Gamma.col(k))/NN;
  }
  return(outPI);
}
arma::mat gmm_standard_mean(arma::mat X, arma::mat Gamma){
  // parameters
  int N = Gamma.n_rows;
  int K = Gamma.n_cols;
  int P = X.n_cols;
  
  // iterate
  arma::mat output(K,P,fill::zeros);
  arma::rowvec tmpvec(P,fill::zeros);
  for (int k=0;k<K;k++){
    tmpvec.fill(0.0);
    for (int n=0;n<N;n++){
      tmpvec += Gamma(n,k)*X.row(n);
    }
    output.row(k) = tmpvec/arma::accu(Gamma.col(k));
  }
  return(output);
}
arma::cube gmm_standard_cov(arma::mat X, arma::mat Gamma, arma::mat Mu, bool usediag=false){
  // parameters
  int N = X.n_rows;
  int P = X.n_cols;
  int K = Mu.n_rows;
  
  // iterate
  arma::cube output(P,P,K,fill::zeros);
  arma::mat  mslice(P,P,fill::zeros);
  arma::rowvec xdiff(P,fill::zeros);
  arma::mat  tmpmat(P,P,fill::zeros);
  double Nk = 0.0;
  for (int k=0;k<K;k++){
    // denominator
    Nk = arma::accu(Gamma.col(k));
    // numerator
    mslice.fill(0.0);
    for (int n=0;n<N;n++){
      xdiff   = X.row(n) - Mu.row(k);
      mslice += Gamma(n,k)*(xdiff.t()*xdiff);
    }
    tmpmat = mslice/Nk;
    if (usediag==true){
      output.slice(k) = arma::diagmat(tmpmat);
    } else {
      output.slice(k) = tmpmat;
    }
  }
  return(output);
}
// [[Rcpp::export]]
Rcpp::List gmm_skeleton(arma::mat& X, int k){
  // PARAMETERS
  int n = X.n_rows; 
  int p = X.n_cols;
  int maxiter = 100;
  
  // PREPARATION : DEFINITION
  arma::mat oldGamma(n,k,fill::zeros);  // (NxK) class membership
  arma::mat newGamma(n,k,fill::zeros);
  arma::mat oldMU(k,p,fill::zeros);     // (KxP) row class means
  arma::mat newMU(k,p,fill::zeros);
  arma::cube oldSIG(p,p,k,fill::zeros); // (PxPxK) covariance slices
  arma::cube newSIG(p,p,k,fill::zeros);
  arma::vec oldPI(k,fill::zeros);       // (K) proportion
  arma::vec newPI(k,fill::zeros);
  double oldLKD = 0.0;
  double newLKD = 0.0;
  
  // PREPARATION : INITIALIZATION
  arma::urowvec initlabel = label_gmm(X, k, 5); // small number of simple task
  for (int i=0; i<n; i++){
    oldGamma(i,initlabel(i)) = 1.0;
  }
  oldPI  = gmm_standard_pi(oldGamma);
  oldMU  = gmm_standard_mean(X, oldGamma);
  oldSIG = gmm_standard_cov(X, oldGamma, oldMU, false); // full covariance
  oldLKD = gmm_loglkd(X, oldPI, oldMU, oldSIG);         // use armadillo routines
  
  // ITERATION
  for (int it=0; it<maxiter; it++){
    // {E} update gamma
    newGamma = gmm_standard_gamma(X, oldMU, oldSIG, oldPI);
    // {M} update parameters + compute loglkd
    newPI  = gmm_standard_pi(newGamma);
    newMU  = gmm_standard_mean(X, newGamma);
    newSIG = gmm_standard_cov(X, newGamma, newMU, false);
    newLKD = gmm_loglkd(X, newPI, newMU, newSIG);
    // Breaking Part
    if ((it>0)&&(newLKD <= oldLKD)){
      break;
    } else {
      oldGamma = newGamma;
      oldPI    = newPI;
      oldMU    = newMU;
      oldSIG   = newSIG;
      oldLKD   = newLKD; 
    }
  }
  
  // USE ARMADILLO ROUTINES
  arma::gmm_full model;
  model.reset(p, k);
  model.set_means(arma::trans(oldMU));
  model.set_fcovs(oldSIG);
  model.set_hefts(arma::trans(oldPI));
  
  // RETURN
  Rcpp::List output;
  output["means"]   = oldMU;
  output["covs"]    = oldSIG;
  output["weight"]  = oldPI;
  output["loglkd"]  = model.sum_log_p(arma::trans(X));
  output["cluster"] = arma::trans(model.assign(arma::trans(X), prob_dist));
  return(output);
}

// (01) gmm_armadillo ==========================================================
// [[Rcpp::export]]
Rcpp::List gmm_armadillo(arma::mat& X, int k, int maxiter, bool usediag){
  if (usediag==true){ // use diagonal
    arma::gmm_diag model;
    bool status = model.learn(arma::trans(X), k, maha_dist, random_subset, 10, maxiter, 1e-11, false);
    if (status==false){
      Rcpp::stop("* gmm : Fitting GMM with diagonal covariance failed.");
    } else {
      // successful, let's take out the elements and return
      
      arma::mat diagcovs = model.dcovs;
      int p = diagcovs.n_rows;
      int k = diagcovs.n_cols;
      arma::cube myfcovs(p,p,k,fill::zeros);
      for (int i=0; i<k; i++){
        myfcovs.slice(i) = arma::diagmat(diagcovs.col(i));
      }
      return Rcpp::List::create(Rcpp::Named("means")=arma::trans(model.means),
                                Rcpp::Named("covs")=myfcovs,
                                Rcpp::Named("weight")=model.hefts,
                                Rcpp::Named("loglkd")=model.sum_log_p(arma::trans(X)),
                                Rcpp::Named("cluster")=arma::trans(model.assign(arma::trans(X), prob_dist)));
    }
  } else {
    arma::gmm_full model;
    bool status = model.learn(arma::trans(X), k, maha_dist, random_subset, 10, maxiter, 1e-11, false);
    if (status==false){
      Rcpp::stop("* gmm : Fitting GMM with full covariance failed.");
    } else {
      // successful, let's take out the elements and return
      return Rcpp::List::create(Rcpp::Named("means")=arma::trans(model.means),
                                Rcpp::Named("covs")=model.fcovs,
                                Rcpp::Named("weight")=model.hefts,
                                Rcpp::Named("loglkd")=model.sum_log_p(arma::trans(X)),
                                Rcpp::Named("cluster")=arma::trans(model.assign(arma::trans(X), prob_dist)));
    }
  }
}
// (02) gmm_11R ================================================================
double gmm11R_objective(arma::mat S, arma::mat X, arma::mat Z, double lambda){
  double term1 = arma::trace(S*X);
  double term2 = log(arma::det(X));
  double term3 = lambda*arma::norm(vectorise(Z),1);
  double output = term1-term2+term3;
  return(output);
}
double gmm11R_shrinkage(double a, double kappa){
  double term1=0;
  double term2=0;
  if (a>kappa){    term1 = a-kappa;  }
  if (a<-kappa){   term2 = -a-kappa; }
  double output = term1-term2;
  return(output);
}
arma::mat gmm11R_ADMMprecision(arma::mat S, double lambda){
  // 1. parameters and set up
  const int max_iter  = 1000;
  const double abstol = 1e-6;
  const double reltol = 1e-3;
  const int    n      = S.n_cols;
  
  arma::mat X(n,n,fill::zeros);
  arma::mat X_hat(n,n,fill::zeros);
  arma::mat Z(n,n,fill::zeros);
  arma::mat Zold(n,n,fill::zeros);
  arma::mat U(n,n,fill::zeros);
  
  double rho   = 1.0;
  double alpha = 1.0;
  
  arma::colvec es(n,fill::zeros);
  arma::colvec xi(n,fill::zeros);
  arma::mat Q(n,n,fill::zeros);
  
  arma::vec objval(max_iter,fill::zeros);
  arma::vec r_norm(max_iter,fill::zeros);
  arma::vec s_norm(max_iter,fill::zeros);
  arma::vec eps_pri(max_iter,fill::zeros);
  arma::vec eps_dual(max_iter,fill::zeros);
  
  for (int k=0;k<max_iter;k++){
    // update X
    eig_sym(es,Q,rho*(Z-U)-S);
    for (int i=0;i<n;i++){
      xi(i) = (es(i)+sqrt(pow(es(i),2)+4*rho))/(2*rho);
    }
    X = Q*arma::diagmat(xi)*Q.t();
    
    // update Z with relaxation
    Zold = Z;
    X_hat = alpha*X + (1-alpha)*Zold;
    for (int i=0;i<n;i++){
      for (int j=0;j<n;j++){
        Z(i,j) = gmm11R_shrinkage(X_hat(i,j)+U(i,j), lambda/rho);
      }
    }
    
    // update U
    U = U + (X_hat-Z);
    
    // diagnostics
    objval(k) = gmm11R_objective(S, X, Z, lambda);
    r_norm(k) = arma::norm(X-Z,"fro");
    s_norm(k) = arma::norm(-rho*(Z-Zold),"fro");
    if (norm(X,"fro")>norm(Z,"fro")){
      eps_pri(k) = static_cast<double>(n)*abstol + reltol*norm(X,"fro");
    } else {
      eps_pri(k) = static_cast<double>(n)*abstol + reltol*norm(Z,"fro");
    }
    eps_dual(k) = static_cast<double>(n)*abstol + reltol*norm(rho*U, "fro");
    
    if ((r_norm(k)<eps_pri(k))&&(s_norm(k)<eps_dual(k))){
      break;
    }
  }
  return(X);
}
arma::cube gmm11R_precision(arma::mat X, arma::mat Gamma, arma::mat Mu, double lambda){
  // parameters
  int N = X.n_rows;
  int P = X.n_cols;
  int K = Mu.n_rows;
  
  arma::mat  Ak(P,P,fill::zeros);
  arma::rowvec xdiff(P,fill::zeros);
  arma::rowvec muk(P,fill::zeros);
  double Nk = 0.0;
  arma::cube output(P,P,K,fill::zeros);
  for (int k=0;k<K;k++){
    // initialize
    Ak.fill(0.0);
    Nk  = arma::accu(Gamma.col(k));
    muk = Mu.row(k);
    for (int n=0;n<N;n++){
      xdiff = X.row(n) - muk;
      Ak   += (Gamma(n,k)/Nk)*(xdiff.t()*xdiff);
    }
    // compute
    output.slice(k) = gmm11R_ADMMprecision(Ak, lambda); // precision to covariance
  }
  return(output);
}
// [[Rcpp::export]]
Rcpp::List gmm_11R(arma::mat& X, int K, double lambda, int maxiter, bool usediag){
  // parameters
  int N = X.n_rows; 
  double NN = static_cast<double>(N);
  int P = X.n_cols;
  int mit = maxiter; // max iter
  
  // preparation : definition
  arma::mat oldGamma(N,K,fill::zeros);  // (NxK) class membership
  arma::mat newGamma(N,K,fill::zeros);
  arma::mat oldMU(K,P,fill::zeros);     // (KxP) row class means
  arma::mat newMU(K,P,fill::zeros);
  arma::cube oldSIG(P,P,K,fill::zeros); // (PxPxK) covariance slices
  arma::cube newSIG(P,P,K,fill::zeros);
  arma::cube oldPRE(P,P,K,fill::zeros);
  arma::cube newPRE(P,P,K,fill::zeros);
  arma::vec oldPI(K,fill::zeros);       // (K) proportion
  arma::vec newPI(K,fill::zeros);
  double oldLKD = 0.0;
  double newLKD = 0.0;
  
  int par0 = (P*K) + (K-1); // mean + proportion
  int parS = 0;             // covariance part
  int oldPar = 0;
  int newPar = 0;
  
  double valijk = 0.0;
  double thres0 = arma::datum::eps;
  
  // preparation : initialize
  arma::urowvec initlabel = label_gmm(X, K, 5); // small number of simple task
  for (int n=0; n<N; n++){
    oldGamma(n,initlabel(n)) = 1.0;
  }
  oldPI  = gmm_standard_pi(oldGamma);
  oldMU  = gmm_standard_mean(X, oldGamma);
  oldPRE = gmm11R_precision(X, oldGamma, oldMU, lambda);
  parS   = 0;
  for (int k=0;k<K;k++){
    for (int i=0;i<P;i++){
      for (int j=i;j<P;j++){
        valijk = std::abs(oldPRE(i,j,k));
        if (valijk >= thres0){
          parS += 1;
        }
      }
    }
  }
  oldPar = par0 + parS;
  for (int k=0;k<K;k++){
    oldSIG.slice(k) = arma::pinv(oldPRE.slice(k));
  }
  if (usediag==true){
    for (int k=0; k<K; k++){
      oldSIG.slice(k) = arma::diagmat(oldSIG.slice(k));
    }
  }
  oldLKD = gmm_loglkd(X, oldPI, oldMU, oldSIG);
  
  // main iteration
  for (int it=0;it<mit;it++){
    // E-step. update Gamma
    newGamma = gmm_standard_gamma(X, oldMU, oldSIG, oldPI);
    // M-step. update parameters + log-likelihood
    newPI  = gmm_standard_pi(newGamma);
    newMU  = gmm_standard_mean(X, newGamma);
    newPRE = gmm11R_precision(X, newGamma, newMU, lambda);
    parS   = 0;
    for (int k=0;k<K;k++){
      for (int i=0;i<P;i++){
        for (int j=i;j<P;j++){
          valijk = std::abs(newPRE(i,j,k));
          if (valijk >= thres0){
            parS += 1;
          }
        }
      }
    }
    newPar = par0 + parS;
    for (int k=0;k<K;k++){
      newSIG.slice(k) = arma::pinv(newPRE.slice(k));
    }
    if (usediag==true){
      for (int k=0; k<K; k++){
        newSIG.slice(k) = arma::diagmat(newSIG.slice(k));
      }
    }
    newLKD = gmm_loglkd(X, newPI, newMU, newSIG);
    // break part
    if ((it>=1)&&(newLKD <= oldLKD)){
      break;
    } else {
      oldGamma = newGamma;
      oldPI    = newPI;
      oldMU    = newMU;
      oldSIG   = newSIG;
      oldLKD   = newLKD;
      oldPRE   = newPRE;
      oldPar   = newPar;
    }
  }
  
  // USE ARMADILLO ROUTINES
  arma::gmm_full model;
  model.reset(P, K);
  model.set_means(arma::trans(oldMU));
  model.set_fcovs(oldSIG);
  model.set_hefts(arma::trans(oldPI));
  
  // RETURN
  Rcpp::List output;
  output["means"]   = oldMU;
  output["covs"]    = oldSIG;
  output["weight"]  = oldPI;
  output["loglkd"]  = model.sum_log_p(arma::trans(X));
  output["cluster"] = arma::trans(model.assign(arma::trans(X), prob_dist));
  return(output);
}

// (03) gmm_16Gfix =============================================================
arma::mat gmm_16Gfix_mean(arma::mat X, arma::mat Gamma, arma::vec Weight){
  // parameters
  int N = Gamma.n_rows;
  int K = Gamma.n_cols;
  int P = X.n_cols;
  
  // iterate
  arma::mat output(K,P,fill::zeros);
  arma::rowvec tmpvec(P,fill::zeros);
  double tmpval=0.0;
  for (int k=0;k<K;k++){
    tmpvec.fill(0.0);
    tmpval = 0.0;
    for (int n=0;n<N;n++){
      tmpvec += Gamma(n,k)*X.row(n)*Weight(n);
      tmpval += Gamma(n,k)*Weight(n);
    }
    output.row(k) = tmpvec/tmpval;
  }
  return(output);
}
arma::cube gmm_16Gfix_cov(arma::mat X, arma::mat Gamma, arma::mat Mu, arma::vec Weight, bool usediag=false){
  // parameters
  int N = X.n_rows;
  int P = X.n_cols;
  int K = Mu.n_rows;
  
  // iterate
  arma::cube output(P,P,K,fill::zeros);
  arma::mat  mslice(P,P,fill::zeros);
  arma::rowvec xdiff(P,fill::zeros);
  arma::mat  tmpmat(P,P,fill::zeros);
  double Nk = 0.0;
  for (int k=0;k<K;k++){
    // denominator
    Nk = arma::accu(Gamma.col(k));
    // numerator
    mslice.fill(0.0);
    for (int n=0;n<N;n++){
      xdiff   = X.row(n) - Mu.row(k);
      mslice += Gamma(n,k)*(xdiff.t()*xdiff)*Weight(n);
    }
    tmpmat = mslice/Nk;
    if (usediag==true){
      output.slice(k) = arma::diagmat(tmpmat);
    } else {
      output.slice(k) = tmpmat;
    }
  }
  return(output);
}
double gmm_16Gfix_loglkd(arma::mat X, arma::vec kPi, arma::mat kMu, arma::cube kSig, arma::vec weight){
  // Parameters
  int n = X.n_rows;
  //int p = X.n_cols;
  int k = kMu.n_rows;
  
  // Compute Elementary
  arma::mat piNN(n,k,fill::zeros);
  for (int i=0; i<n; i++){
    for (int j=0; j<k; j++){
      piNN(i,j) = kPi(j)*eval_gaussian_single(X.row(i), kMu.row(j), kSig.slice(j)/weight(i), false);
    }
  }
  
  // Compute
  arma::vec before_log(n,fill::zeros);
  for (int i=0; i<n; i++){
    before_log(i) = std::log(arma::accu(piNN.row(i)));
  }
  
  // Return
  return(arma::accu(before_log));
}
arma::uvec gmm_16Gfix_label(arma::mat& X, arma::mat parMU, arma::cube parSIG, arma::vec parPI, arma::vec weight){
  // parameters
  int N = X.n_rows;
  int K = parSIG.n_slices;
  // compute gamma
  arma::mat piNN(N,K,fill::zeros);
  for (int n=0; n<N; n++){
    for (int k=0; k<K; k++){
      piNN(n,k) = parPI(k)*eval_gaussian_single(X.row(n), parMU.row(k), parSIG.slice(k)/weight(n), false);
    }
  }
  // for each component, find the maximal
  arma::uvec output(N,fill::zeros);
  for (int n=0; n<N; n++){
    output(n) = arma::index_max(piNN.row(n));
  }
  return(output);
}

// [[Rcpp::export]]
Rcpp::List gmm_16Gfix(arma::mat& X, int k, arma::vec weight, int maxiter, bool usediag){
  // PARAMETERS
  int n = X.n_rows; 
  int p = X.n_cols;
  
  // PREPARATION : DEFINITION
  arma::mat oldGamma(n,k,fill::zeros);  // (NxK) class membership
  arma::mat newGamma(n,k,fill::zeros);
  arma::mat oldMU(k,p,fill::zeros);     // (KxP) row class means
  arma::mat newMU(k,p,fill::zeros);
  arma::cube oldSIG(p,p,k,fill::zeros); // (PxPxK) covariance slices
  arma::cube newSIG(p,p,k,fill::zeros);
  arma::vec oldPI(k,fill::zeros);       // (K) proportion
  arma::vec newPI(k,fill::zeros);
  double oldLKD = 0.0;
  double newLKD = 0.0;
  
  // PREPARATION : INITIALIZATION
  arma::urowvec initlabel = label_gmm(X, k, 5); // small number of simple task
  for (int i=0; i<n; i++){
    oldGamma(i,initlabel(i)) = 1.0;
  }
  oldPI  = gmm_standard_pi(oldGamma);
  oldMU  = gmm_16Gfix_mean(X, oldGamma, weight);
  oldSIG = gmm_16Gfix_cov(X, oldGamma, oldMU, weight, usediag);
  oldLKD = gmm_16Gfix_loglkd(X, oldPI, oldMU, oldSIG, weight);
  
  // ITERATION
  for (int it=0; it<maxiter; it++){
    // {E} update gamma
    newGamma = gmm_standard_gamma(X, oldMU, oldSIG, oldPI);
    // {M} update parameters + compute loglkd
    newPI  = gmm_standard_pi(newGamma);
    newMU  = gmm_16Gfix_mean(X, newGamma, weight);
    newSIG = gmm_16Gfix_cov(X, newGamma, newMU, weight, usediag);
    newLKD = gmm_16Gfix_loglkd(X, newPI, newMU, newSIG, weight);
    // Breaking Part
    if ((it>0)&&(newLKD <= oldLKD)){
      break;
    } else {
      oldGamma = newGamma;
      oldPI    = newPI;
      oldMU    = newMU;
      oldSIG   = newSIG;
      oldLKD   = newLKD; 
    }
  }
  
  // USE ARMADILLO ROUTINES
  arma::gmm_full model;
  model.reset(p, k);
  model.set_means(arma::trans(oldMU));
  model.set_fcovs(oldSIG);
  model.set_hefts(arma::trans(oldPI));
  
  // RETURN
  Rcpp::List output;
  output["means"]   = oldMU;
  output["covs"]    = oldSIG;
  output["weight"]  = oldPI;
  output["loglkd"]  = model.sum_log_p(arma::trans(X));
  output["cluster"] = gmm_16Gfix_label(X, oldMU, oldSIG, oldPI, weight);
  return(output);
}