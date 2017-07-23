mult_tnorm_surv <- function(x, mean, sd, limits, twosided = TRUE)
{
  
  limits <- limits[order(limits[,1]),]
  pnormstd <- function(x) pnorm((x - mean)/sd)
  limitX <- which(limits[,2] > x & limits[,1] < x)
  
  # allow for multiple means
  if(length(mean)==1)
  {
    
    denom <- 
      sum(sapply(1:nrow(limits), function(i) 
        pnormstd(limits[i,2]) - pnormstd(limits[i,1])))
    
  }else{
    
    denom <- sapply(1:nrow(limits), function(i) 
      pnormstd(limits[i,2]) - pnormstd(limits[i,1]))
    denom <- apply(denom, 1, sum)
    
  }
  
  limitXb <- limits[limitX,2]
  
  if(length(mean)==1)
  {
    
    massLargerXb <- if(nrow(limits)>limitX) 
      sum(sapply((limitX+1):nrow(limits), function(i) 
        pnormstd(limits[i,2]) - pnormstd(limits[i,1]))) else 0
    
  }else{
    
    if(nrow(limits)>limitX)
    {
      
      massLargerXb <- sapply((limitX+1):nrow(limits), function(i) 
        pnormstd(limits[i,2]) - pnormstd(limits[i,1]))
      massLargerXb <- apply(massLargerXb, 1, sum)
      
    }else{
      
      massLargerXb <- rep(0, length(mean))
      
    }
    
  }
  
  nom <- (massLargerXb + pnormstd(limitXb)) - pnormstd(x)
  
  this_pval <- nom / denom
  this_pval[is.nan(this_pval) & mean < x] <- 0
  this_pval[is.nan(this_pval) & mean > x] <- 1
  
  if(all(is.na(this_pval) | is.nan(this_pval)))
    this_pval <- bryc.tnorm.surv(z = x, sd = sd, 
                                 a = limits[limitX,1], 
                                 b = limits[limitX,2])
  
  if(twosided) this_pval <- sapply(this_pval, 
                                   function(x) 2*min(c(x, 1-x)))
  
  return(this_pval)
  
}


#' Calculate confidence interval for truncated normal variable
#' 
#' @param meanEst mean
#' @param sd standard deviation
#' @param limits \code{Intervals} object containing all limits
#' @param alpha value specifying the interval coverage, which is \code{(1-alpha)}
#' @param gridpts,griddepth parameters for grid search used in the confidence interval
#' function of the package \code{selectiveInference}.
#' @param range gives the range \code{range*sd}, in which the quantile values are searched
#' 
#' @examples 
#' library(intervals)
#' # one restriction -- in this case equal to a call to \code{\link{msm::qtnorm}}
#' int1 <- Intervals(c(-3,1))
#' mean = -0.2
#' sd = 1
#' ci_tnorm(meanEst = mean, sd = sd, limits = int1, alpha = 0.05)
#' 
#' # two intervals -- 1-alpha coverage not possible
#' (int2 <- Intervals(matrix(c(-3,-2,-1,1), 
#'                           ncol = 2, byrow = TRUE)))
#' ci_tnorm(meanEst = mean, sd = sd, limits = int2, alpha = 0.05)
#' 
#' # two intervals -- 1-alpha coverage possible, but interval not symmetric
#' (int3 <- Intervals(matrix(c(-3,-2.5,-1.9,1), 
#'                           ncol = 2, byrow =  TRUE)))
#' ci_tnorm(meanEst = mean, sd = sd, limits = int3, alpha = 0.05)
#' 
#' @importFrom msm qtnorm
#' @export
#' 
ci_tnorm <- function(meanEst, sd, limits, alpha, 
                     gridpts = 500, griddepth = 3, range = 100)
{
  
  # grid search adapted from the selectiveInference package
  
  limitsX <- which(limits[,2] > meanEst & limits[,1] < meanEst)
  
  xg = seq(-range * sd, range * sd, length = gridpts)
  
  if(nrow(limits)==1){
    
    fun = function(x) { 
      tnorm.surv(z = meanEst, mean = x, sd = sd, 
                 a = limits[limitsX, 1],
                 b = limits[limitsX, 2]) }
    
  }else{
    
    fun = function(x) { 
      mult_tnorm_surv(x = meanEst, mean = x, sd = sd, 
                      limits = limits, twosided = FALSE) }
    
  }
  
  int = grid.search(xg, fun, alpha/2, 1-alpha/2, 
                    gridpts, griddepth)
  
  return(Intervals(int))
  
  
}
