#' Median of GMMs under \eqn{L^2} Metric
#' 
#' 
#' 
#' 
#' @examples 
#' # -------------------------------------------------------------
#' #                   Mixtures for SMILEY Data
#' #
#' # STEP 1. Generate 10 datasets with noise
#' # STEP 2. Fit GMM with varying number of components (k=3-8)
#' # STEP 3. Combine GMMs with uniform weights and Visualize
#' # -------------------------------------------------------------
#' # STEP 1. 10 datasets with noise
#' list_data = list()
#' for (i in 1:10){
#'   list_data[[i]] = T4cluster::genSMILEY(sd=0.25)
#' }
#' 
#' # STEP 2. Fit GMM with varying k
#' list_gmm = list()
#' for (i in 1:10){
#'   list_gmm[[i]] = gmmtools::gmm(list_data[[i]]$data, k=sample(3:8, 1))
#' }
#' 
#' # STEP 3. Find Median of GMMs under L2 metric
#' medgmm = medianL2(list_gmm)
#' 
#' # prepare grid and density evaluation
#' grid.data = expand.grid(x=seq(from=-1.5,to=1.5,length.out=250),
#'                         y=seq(from=-1.5,to=1.5,length.out=250))
#' grid.prob = gmmdensity(medgmm, grid.data)
#' grid.df   = as.data.frame(cbind(grid.data, density=grid.prob))
#' 
#' \dontrun{
#' # we use 'ggplot2' for visualization
#' require(ggplot2)
#' ggplot(grid.df, aes(x=x,y=y,z=density)) + 
#'   geom_raster(aes(fill=density)) + 
#'   geom_contour(colour="white") +
#'   scale_fill_viridis_c() +
#'   scale_x_continuous(expand = c(0, 0)) +
#'   scale_y_continuous(expand = c(0, 0)) +
#'   coord_fixed(xlim=c(-1.5,1.5), ylim=c(-1.5,1.5)) + 
#'   ggtitle("L2 Median of 10 GMMs")
#' }
#' 
#' @concept combine
#' @export
medianL2 <- function(gmmlist, weight=NULL, ...){
  ## INPUTS : EXPLICIT
  check_gmmlist(gmmlist, "medianL2")
  K = length(gmmlist)
  if ((length(weight)<1)&&(is.null(NULL))){
    weight = rep(1/K, K)
  } else {
    cond1 = is.vector(weight)
    cond2 = (length(weight)==K)
    cond3 = all(weight > 0)
    if (!(cond1&&cond2&&cond3)){
      stop(paste0("* medianL2 : input 'weight' should be a vector of length ",K," with nonnegative weights."))
    }
    weight = weight/base::sum(weight)
  }
  
  ## INPUTS : IMPLICIT
  mypars  = list(...)
  pnames  = names(mypars)
  if ("maxiter"%in%pnames){
    maxiter = max(round(mypars$maxiter), 5)
  } else {
    maxiter = 50
  }
  if ("abstol"%in%pnames){
    abstol = max(sqrt(.Machine$double.eps), as.double(mypars$abstol))
  } else {
    abstol = 1e-6
  }
  
  ## RUN C++ ROUTINE
  cpprun  = cpp_combine_medianL2(gmmlist, weight, maxiter, abstol)
  
  ## PREPARE AND RETURN
  output = list()
  output$mean     = cpprun$means
  output$variance = cpprun$covs
  output$weight   = as.vector(cpprun$weight)
  return(structure(output, class="gmmtools"))
}