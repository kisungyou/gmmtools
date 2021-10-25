#' \eqn{L_2} Distance
#' 
#' \eqn{L_2} distance between two measures \eqn{P} and \eqn{Q} with corresponding 
#' probability density functions \eqn{p(x)} and \eqn{q(x)} is defined as
#' \deqn{
#' D_{L_2}(P,Q) = \sqrt{\int \left(p(x)-q(x)\right)^2 dx}
#' }
#' which admits closed-form evaluation when \eqn{p(x)} (\code{gmmobj1}) and \eqn{q(x)} (\code{gmmobj2}) are two mixtures of Gaussians. 
#' 
#' @param gmm1 a \code{gmmtools} object.
#' @param gmm2 a \code{gmmtools} object.
#' 
#' @return computed distance value.
#' 
#' @examples 
#' \donttest{
#' # -------------------------------------------------------------
#' #                   Distance for Gaussian Mixtures
#' #
#' # Data 1 : use SMILEY data 'gensmiley()' function.
#' # Data 2 : SMILEY data is translated +3 and rotated.
#' # Data 3 : SMILEY data is translated +9 and rotated.
#' # -------------------------------------------------------------
#' ## GENERATE DATA
#' #  set up
#' ndata = 10
#' ntot  = 3*ndata
#' rot   = qr.Q(qr(matrix(rnorm(4),ncol=2)))
#' 
#' #  generate
#' list_data = list()
#' for (i in 1:ndata){
#'   list_data[[i]]           = (T4cluster::genSMILEY(n=150, sd=0.1)$data)
#'   list_data[[i+ndata]]     = (T4cluster::genSMILEY(n=150, sd=0.1)$data)%*%rot + 3
#'   list_data[[i+(2*ndata)]] = (T4cluster::genSMILEY(n=150, sd=0.1)$data)%*%rot + 9
#' }
#' 
#' ## FIT GMM MODELS WITH K=5 : Just Arbitrary Choice
#' list_gmm = list()
#' for (i in 1:ntot){
#'   list_gmm[[i]] = gmm(list_data[[i]], k=5)
#' }
#' 
#' ## COMPUTE PAIRWISE DISTANCE
#' distL2 = array(0,c(ntot,ntot))
#' for (i in 1:(ntot-1)){
#'   gi = list_gmm[[i]]
#'   for (j in (i+1):ntot){
#'     gj = list_gmm[[j]]
#'     distL2[i,j] <- distL2[j,i] <- dissL2(gi, gj)
#'   }
#' }
#' 
#' ## VISUALIZE
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2), pty="s")
#' plot(list_data[[1]], main="SMILEY data", pch=19)
#' image(distL2[,ntot:1], axes=FALSE, col=gray((0:64)/64),
#'       main="L2 Distance")
#' par(opar)
#' }
#' 
#' @concept dissimilarity
#' @export
dissL2 <- function(gmm1, gmm2){
  ## INPUTS
  check_gmmobj(gmm1, "dissL2")
  check_gmmobj(gmm2, "dissL2") 

  ## COMPUTE WITH C++
  output = as.double(cpp_gmmdist_l2(gmm1$weight, gmm1$mean, gmm1$variance, 
                                    gmm2$weight, gmm2$mean, gmm2$variance))
  return(output)
}

# # generate
# nobj = 5
# displace = 10
# list_data = list()
# for (i in 1:nobj){
#   list_data[[i]]    = (T4cluster::genSMILEY(n=150, sd=0.1)$data)
#   list_data[[i+nobj]] = (T4cluster::genSMILEY(n=150, sd=0.1)$data) + displace
# }
# list_gmm = list()
# for (i in 1:(2*nobj)){
#   list_gmm[[i]] = gmmtools::gmm(list_data[[i]], k=5)
# }
# dmat = array(0,c(2*nobj,2*nobj))
# for (i in 1:(2*nobj - 1)){
#   for (j in (i+1):(2*nobj)){
#     dmat[i,j] <- dmat[j,i] <- dissL2(list_gmm[[i]], list_gmm[[j]])
#   }
# }
# image(dmat)


# tmpdat = (T4cluster::genSMILEY(n=150, sd=0.1)$data)
# hey1   = gmm(tmpdat, k=5)
# hey2   = gmm(tmpdat+5, k=5)
# 
# xgrid = seq(from=0.0001, to=100, length.out=50)
# ygrid = rep(0, 50)
# for (i in 1:50){
#   ygrid[i] = dissHD(hey1, hey2, theta=xgrid[i])
#   print(paste0("* iteration ",i, " complete."))
# }
# plot(xgrid, ygrid, "b")
# abline(h=dissL2(hey1, hey2), col="red")

# single_gaussian <- function(x, mu, sigma){
#   d = length(mu)
#   xmu = x-mu
#   
#   term1 = -sum(as.vector(solve(sigma,xmu))*xmu)
#   term2 = -(d/2)*log(2*pi)
#   term3 = -0.5*log(det(sigma))
#   return(exp(term1+term2+term3))
# }
# 
# p = 3
# par_mu  = rnorm(p)
# par_sig = cov(matrix(rnorm(10*p),ncol=p))
# 
# vec_theta = 10^(seq(from=-4,to=1, length.out=100))
# vec_vals  = rep(0,100)
# for (i in 1:100){
#   now_sig = (2*par_sig) + diag(p)*vec_theta[i]
#   vec_vals[i] = single_gaussian(par_mu, par_mu, now_sig)*((vec_theta[i])^p)
# }