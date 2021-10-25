#' Wrap Multivariate Gaussian Mixture Model
#' 
#' Given a set of user-specified parameters, this function returns an object of 
#' S3 class \code{gmmtools} that represents an \eqn{n}-dimensional Gaussian mixture model. 
#' 
#' @param weight a length-\eqn{K} vector of class weights that sum to 1.
#' @param mean a \eqn{(K\times n)} matrix of component means.
#' @param variance a \eqn{(n\times n\times K)} array of component covariances.
#' 
#' @return a S3 object of \code{gmmtools} class.
#' 
#' @examples 
#' \donttest{
#' # a simple example with 2 components in R^3
#' my_pi  = c(0.4, 0.6)
#' my_mu  = matrix(rnorm(2*3), nrow=2)
#' my_var = array(0,c(3,3,2))
#' for (i in 1:2){
#'   tmpmat = matrix(rnorm(10*3),ncol=3)
#'   my_var[,,i] = stats::cov(tmpmat)
#' }
#' 
#' # wrap with given parameters
#' gmm3d  = wrapmd(my_pi, my_mu, my_var)
#' }
#' 
#' @concept util
#' @export
wrapmd <- function(weight, mean, variance){
  ## Inputs
  myweight = as.vector(weight)
  myweight = myweight/base::sum(myweight)
  if (any(myweight<=0)){
    stop("* wrapmd : 'weight' must be a vector of positive real numbers that sum to 1.")
  }
  myk = length(myweight)
  
  ## Case Branching
  if (myk < 2){
    tmpmean = as.vector(mean)
    tmpvars = as.matrix(variance)
    p = length(tmpmean)
    
    mymean = matrix(tmpmean, nrow=1)
    myvars = array(0,c(p,p,1))
    myvars[,,1] = tmpvars
  } else {
    # take care of means
    tmpmean = as.matrix(mean)
    if (nrow(tmpmean)!=myk){
      stop(paste0("* wrapmd : 'mean' should be a matrix of ",myk," rows."))
    }
    p = ncol(tmpmean)
    mymean = tmpmean
    
    # and variances
    myvars = array(0,c(p,p,myk))
    if (length(dim(variance))!=3){
      stop(paste0("* wrapmd : 'variance' should be a 3d-array of ",myk," slices."))
    }
    for (i in 1:myk){
      tgtmat = as.matrix(variance[,,i])
      if (!check_sympd(tgtmat)){
        stop(paste0("* wrapmd : ",i,"-th covariance matrix is not symmetric and positive-definite."))
      } else {
        myvars[,,i] = variance[,,i]  
      }
    }
  }
  
  ## WRAP AND RETURN
  output = list()
  output$mean     = mymean
  output$variance = myvars
  output$weight   = myweight
  return(structure(output, class="gmmtools"))  
}