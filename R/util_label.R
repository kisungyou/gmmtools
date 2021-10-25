#' Label Prediction of Data given a Fitted GMM Model
#' 
#' Given new data with the fitted GMM models, predict the cluster labels for the newly provided data according to the fitted model. Please note that this function benefits from the coherent structures of the package itself so that the input \code{gmmobj} must be an output of any GMM-type routines in \pkg{gmmtools} package.
#' 
#' @param gmm an output of any GMM routines in our package of \code{gmmtools} class.
#' @param data an \eqn{(n\times p)} query points.
#' 
#' @return a length-\eqn{n} vector of class labels.
#' 
#' @examples
#' # -------------------------------------------------------------
#' #            clustering with 'iris' dataset
#' # -------------------------------------------------------------
#' ## PREPARE
#' data(iris)
#' X   = as.matrix(iris[,1:4])
#' 
#' ## LEARN WITH TESTING DATA FOR K=3
#' cl3 = gmm(X, k=3)
#' 
#' ## PREDICT LABEL FOR THE SAME DATA 
#' lab.test = gmmlabel(cl3, X)
#' lab.pred = cl3$cluster
#' 
#' ## EMBEDDING WITH PCA
#' X2d    = Rdimtools::do.pca(X, ndim=2)$Y
#' 
#' ## VISUALIZATION
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2), pty="s")
#' plot(X2d, col=lab.pred, pch=19, main="original results")
#' plot(X2d, col=lab.test, pch=19, main="predicted label")
#' par(opar)
#' 
#' @concept util
#' @export
gmmlabel <- function(gmm, data){
  # PREPARE
  check_gmmobj(gmm, "gmmdensity")
  mydata = check_inputmat(data)
  
  # RUN
  output = eval_label(mydata, gmm$mean, gmm$variance, gmm$weight)
  return(as.vector(output)+1)
}
