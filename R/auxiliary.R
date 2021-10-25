## AUXILIARY FUNCTIONS
#  check_inputmat : return output as row-stacked matrix
#  check_sympd    : check whether a matrix is symmetric and positive-definite
#  check_gmmobj   : whether an input is a valid GMM object
#  check_gmmlist  : whether a list of GMM objects




# check_inputmat ----------------------------------------------------------
#' @keywords internal
#' @noRd
check_inputmat <- function(x){
  if (is.vector(x)){
    output = matrix(x, ncol=1)
    return(output)
  } else {
    if (is.matrix(x)){
      return(x)
    } else {
      stop("* input should be either a vector or a matrix.")  
    }
  }
}


# check_sympd -------------------------------------------------------------
#' @keywords internal
#' @noRd
check_sympd <- function(tmpmat){
  cond1 = isSymmetric(tmpmat)
  cond2 = (min(base::eigen(tmpmat)$values) > .Machine$double.eps)
  if (cond1&&cond2){
    return(TRUE)
  } else {
    return(FALSE)
  }
}


# check_gmmobj ------------------------------------------------------------
#' @keywords internal
#' @noRd
check_gmmobj <- function(gmmobj, fname){
  if (!inherits(gmmobj,"gmmtools")){
    stop(paste0("* ",fname," : input 'gmm' should be an object from any gmm functions in our package."))
  }
  objnames = names(gmmobj)
  cond1 = ("mean"%in%objnames)
  cond2 = ("variance"%in%objnames)
  cond3 = ("weight"%in%objnames)
  if (!(cond1&&cond2&&cond3)){
    stop(paste0("* ",fname," : input 'gmm' is not a valid object with mean, variance, and weight parameters."))
  } 
}

# check_gmmlist -----------------------------------------------------------
#' @keywords internal
#' @noRd
check_gmmlist <- function(gmmlist, fname){
  if (!is.list(gmmlist)){
    stop(paste0("* ",fname," : input 'gmmlist' is not a list."))
  } else {
    K    = length(gmmlist)
    veck = rep(0,K)
    for (k in 1:K){
      tgtgmm = gmmlist[[k]]
      objnames = names(tgtgmm)
      
      cond1 = inherits(tgtgmm,"gmmtools")
      cond2 = ("mean"%in%objnames)
      cond3 = ("variance"%in%objnames)
      cond4 = ("weight"%in%objnames)
      
      if (!(cond1&&cond2&&cond3&&cond4)){
        stop(paste0("* ",fname," : input 'gmmobj' contains a non-GMM-like object."))
      }
      veck[k] = base::ncol(tgtgmm$mean)
    }
    if (length(unique(veck))!=1){
      stop(paste0("* ",fname," : input 'gmmobj' consists of models fitted in different dimension."))
    }
  }
}

