#' Draw Samples from a Fitted GMM Model
#' 
#' Given a fitted GMM model in \eqn{\mathbf{R}^p}, it generates \eqn{n} samples 
#' according to the fitted model. 
#' 
#' @param gmm an output of any GMM routines in our package of \code{gmmtools} class.
#' @param n the number of observations to be drawn.
#' 
#' @return an \eqn{(n\times p)} generated data matrix.
#' 
#' @examples
#' # -------------------------------------------------------------
#' #               example using 'iris' dataset
#' # -------------------------------------------------------------
#' ## PREPARE
#' data(iris)
#' X   = as.matrix(iris[,1:4])
#' 
#' ## FIT THE MODEL WITH K=4
#' cl3 = gmm(X, k=4)
#' 
#' ## GENERATE 150 SAMPLES FROM THE MODEL
#' Z = gmmsample(cl3, 150)
#' 
#' ## PCA FOR ORIGINAL AND GENERATED DATA
#' X2trf = Rdimtools::do.pca(X, ndim=2)
#' X2ori = X2trf$Y
#' X2gen = (Z-matrix(rep(colMeans(Z),150),ncol=4,byrow=TRUE))%*%X2trf$projection
#' 
#' ## VISUALIZE
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2), pty="s")
#' plot(X2ori, pch=19, cex=0.8, main="original data")
#' plot(X2gen, pch=19, cex=0.8, main="generated data")
#' par(opar)
#' 
#' @concept util
#' @export
gmmsample <- function(gmm, n){
  ## PREPARE : EXPLICIT INPUT
  myn = max(1, round(n))
  check_gmmobj(gmm,"gmmsample")
  
  ## RUN
  output = gmm_sample(myn, gmm$weight, gmm$mean, gmm$variance)
  return(output)
}