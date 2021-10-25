#' Wrap Univariate Gaussian Mixture Model
#' 
#' Given a set of user-specified parameters, this function returns an object of 
#' S3 class \code{gmmtools} that represents a 1-dimensional Gaussian mixture model. 
#' 
#' @param weight a length-\eqn{K} vector of class weights that sum to 1.
#' @param mean a length-\eqn{K} vector of component means.
#' @param variance a length-\eqn{K} vector of component variances.
#' 
#' @return a S3 object of \code{gmmtools} class.
#' 
#' @examples 
#' \donttest{
#' # a simple example with 2 components
#' my_pi  = c(0.4, 0.6)
#' my_mu  = c(-1, 1)
#' my_var = c(2, 2)
#' 
#' # wrap with given parameters
#' gmm1d  = wrap1d(my_pi, my_mu, my_var)
#' }
#' 
#' @concept util
#' @export
wrap1d <- function(weight, mean, variance){
  ## Inputs
  myweight = as.vector(weight)
  myweight = myweight/base::sum(myweight)
  mymean   = as.vector(mean)
  myvars   = as.vector(variance)
  myk      = length(myweight)
  
  ## Checkers
  if (any(myweight<=0)){
    stop("* wrap1d : 'weight' must be a vector of positive real numbers that sum to 1.")
  }
  if (length(mymean)!=myk){
    stop(paste0("* wrap1d : 'mean' must be a vector of length ",myk,"."))
  }
  if (length(myvars)!=myk){
    stop(paste0("* wrap1d : 'variance' must be a vector of length ",myk,"."))
  }
  if (any(myvars<=0)){
    stop("* wrap1d : 'variance' must be a vector of positive real numbers.")
  }
  
  ## Transform into 'gmmtools' 
  gmm_mean = matrix(mymean)
  gmm_vars = array(0,c(1,1,myk))
  for (k in 1:myk){
    gmm_vars[,,k] = myvars[k]
  }
  
  
  ## WRAP AND RETURN
  output = list()
  output$mean     = gmm_mean
  output$variance = gmm_vars
  output$weight   = myweight
  return(structure(output, class="gmmtools"))  
}