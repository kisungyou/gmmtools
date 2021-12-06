#' GMM Reduction by You (2022)
#' 
#' Given a GMM model with \eqn{N} components for large \eqn{N},
#' model reduction is executed in a graph-theoretic way so that it involves 
#' no tuning parameters. For example, this function does not require a 
#' desired number of components for a compressed mixtures since it is automatically 
#' found from the data.
#' 
#' @param gmm a \code{gmmtools} object of \eqn{N} components.
#' @param merger the type of merging machinery, either \code{"moment"} or \code{"wasserstein"}.
#' 
#' @return a \code{gmmtools} object of \eqn{M} components, where \eqn{M} is 
#' automatically selected via graph clustering algorithm.
#' 
#' @examples 
#' # -------------------------------------------------------------
#' #              Reduction of Mixture for SMILEY Data
#' #
#' # From multiple SMILEY data, we fit 7-component GMM for each 
#' # data and average of models is taken to give a GMM fit with 
#' # large number of redundant components.
#' # -------------------------------------------------------------
#' # Generate 20 datasets with noise and fix GMM 
#' list_gmm  = list()
#' for (i in 1:20){
#'   data_i = T4cluster::genSMILEY(sd=0.25)$data
#'   list_gmm[[i]] = gmm(data_i, k=7)
#' }
#' 
#' # Find the average of models
#' gcenter = wsum(list_gmm)
#' 
#' # Do Reduction using Two Methods
#' grM = gmr2022Y(gcenter, merger="moment")
#' grW = gmr2022Y(gcenter, merger="wasserstein")
#' 
#' \dontrun{
#' # WARNING : RUN THIS PART FOR VISUALIZATION WITH GGPLOT2
#' # prepare grid and density evaluation
#' require("ggplot2")
#' npts  = 200
#' pgrid = as.matrix(expand.grid(x=seq(from=-1.5,to=1.5,length.out=npts),
#'                   y=seq(from=-1.5,to=1.5,length.out=npts)))
#' probM = gmmdensity(grM, data=pgrid)
#' probW = gmmdensity(grW, data=pgrid)
#'   
#' # wrap as a single dataframe
#' obj1 = rbind(cbind(pgrid, probM), cbind(pgrid, probW))
#' obj2 = as.factor(rep(c(1,2),each=(npts^2)))
#' odf  = data.frame(x=obj1[,1], y=obj1[,2], density=obj1[,3], class=obj2)
#' levels(odf$class) = c("moment","wasserstein")
#' 
#' # plot!
#' ggplot2::ggplot(odf, aes(x=x,y=y,z=density)) + 
#'   facet_grid(. ~ class) + 
#'   geom_raster(aes(fill=density)) + 
#'   geom_contour(colour="white") +
#'   scale_fill_viridis_c() +
#'   scale_x_continuous(expand = c(0, 0)) +
#'   scale_y_continuous(expand = c(0, 0)) +
#'   coord_fixed(xlim=c(-1.5,1.5), ylim=c(-1.5,1.5)) + 
#'   ggtitle("Reduction of 140-Component GMM via You (2022).")
#' }  
#' 
#' @concept reduction
#' @export
gmr2022Y <- function(gmm, merger=c("moment","wasserstein")){
  ## INPUTS ====================================================================
  check_gmmobj(gmm, "gmr2022Y")
  merger = match.arg(merger)
  
  ## MAIN COMPUTATION ==========================================================
  #  STEP 1. BHATTACHARYYA DISTANCE COMPUTATION
  mat_BD = interdist_bhat(gmm$mean, gmm$variance)

  #  STEP 2. BUILD AN AFFINITY & IGRAPH OBJECT
  mat_A   = 1-(base::acos(exp(-mat_BD))*2/pi)
  graph_A = igraph::graph_from_adjacency_matrix(mat_A, mode="undirected", weighted=TRUE, diag=FALSE)
  
  #  STEP 3. MODULARITY OPTIMIZATION
  louvain_opt = igraph::cluster_louvain(graph_A)
  louvain_lab = igraph::membership(louvain_opt)
  lab_unique  = unique(louvain_lab)
  
  #  STEP 4. MERGE
  if (length(lab_unique)<2){ # case : singleton
    # merge
    if (all(merger=="moment")){
      cpprun = cpp_collapse_MPM(gmm$weight, gmm$mean, gmm$variance)
    } else {
      cpprun = cpp_collapse_W2(gmm$weight, gmm$mean, gmm$variance)
    }
    # wrap 
    p       = base::ncol(gmm$mean)
    outmean = matrix(as.vector(cpprun$mean), nrow=1)
    outvars = array(0,c(p,p,1))
    outvars[,,1] = as.matrix(cpprun$variance)
    
    output = list()
    output$weight   = 1
    output$mean     = outmean
    output$variance = outvars
    return(structure(output, class="gmmtools"))  
  } else {                   # case : 2 or more clusters
    # setup's
    p = base::ncol(gmm$mean)
    K = length(lab_unique)
    
    output = list()
    output$weight   = rep(0,K)
    output$mean     = array(0,c(K,p))
    output$variance = array(0,c(p,p,K))
    
    # compute & return
    for (k in 1:K){
      idk = which(louvain_lab==lab_unique[k])
      if (length(idk) < 2){
        output$weight[k] = gmm$weight[idk]
        output$mean[k,]  = as.vector(gmm$mean[idk,])
        output$variance[,,k] = as.matrix(gmm$variance[,,idk])
      } else {
        part_weight = gmm$weight[idk]
        part_weight = part_weight/base::sum(part_weight)
        if (all(merger=="moment")){
          cpprun = cpp_collapse_MPM(part_weight, gmm$mean[idk,], gmm$variance[,,idk])
        } else {
          cpprun = cpp_collapse_W2(part_weight, gmm$mean[idk,], gmm$variance[,,idk])
        }
        output$weight[k]     = base::sum(gmm$weight[idk])
        output$mean[k,]      = as.vector(cpprun$mean)
        output$variance[,,k] = as.matrix(cpprun$variance)
      }
    }
    return(structure(output, class="gmmtools"))  
  }
}