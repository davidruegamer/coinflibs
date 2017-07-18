mult_tnorm_surv <- function(x, mean, sd, limits, lower = FALSE, twosided = TRUE)
{
  
  limits <- limits[order(limits[,1]),]
  pnormstd <- function(x) pnorm((x - mean)/sd)
  limitX <- which(limits[,2] > x & limits[,1] < x)
  
  denom <- sum(sapply(1:nrow(limits), function(i) pnormstd(limits[i,2]) - pnormstd(limits[i,1])))
  
  if(lower){
    
    limitXa <- limits[limitX,1]
    massSmallerXa <- if(limitX!=1) 
      sum(sapply(1:(limitX-1), function(i) 
        pnormstd(limits[i,2]) - pnormstd(limits[i,1]))) else 0
    
    nom <- pnormstd(x) - (massSmallerXa + pnormstd(limitXa))
    
  }else{
    
    limitXb <- limits[limitX,2]
    massLargerXb <- if(nrow(limits)>limitX) 
      sum(sapply((limitX+1):nrow(limits), function(i) 
        pnormstd(limits[i,2]) - pnormstd(limits[i,1]))) else 0
    
    nom <- (massLargerXb + pnormstd(limitXb)) - pnormstd(x)
    
  }
  
  this_pval <- nom / denom
  
  if(all(is.na(this_pval) | is.nan(this_pval)) & !lower)
    this_pval <- bryc.tnorm.surv(z = x, sd = sd, 
                                 a = limits[limitX,1], 
                                 b = limits[limitX,2])
  
  if(twosided) this_pval <- 2*min(c(this_pval, 1-this_pval))
  
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
ci_tnorm <- function(meanEst, sd, limits, alpha, gridpts=200, griddepth=3)
{
  
  # limits <- limits[order(limits[,1]),]
  # pnormstd <- function(x) pnorm((x - meanEst)/sd)
  # limitX <- which(limits[,2] > meanEst & limits[,1] < meanEst)
  # 
  # limitXa <- limits[limitX,1]
  # massSmallerXa <- if(limitX!=1) 
  #   sum(sapply(1:(limitX-1), function(i) 
  #     pnormstd(limits[i,2]) - pnormstd(limits[i,1]))) else 0
  # 
  # nom <- (massSmallerXa + pnormstd(limitXa))
  # denom <- sum(sapply(1:nrow(limits), function(i) 
  #   pnormstd(limits[i,2]) - pnormstd(limits[i,1])))
  # 
  # massUpper <- (pnormstd(limits[limitX,2]) - nom) / denom
  # massLower <- (nom - pnormstd(limits[limitX,1])) / denom
  # 
  # thisQfun <- function(p) qtnorm(p, mean = meanH0, sd = sd, 
  #                                lower = limits[limitX, 1],
  #                                upper = limits[limitX, 2])
  # thisCov <- massUpper - massLower
  # 
  # ### case discrimination
  # 
  # if(thisCov > 1-alpha){
  #   
  #   lowdef <- massLower > (alpha/2)
  #   updef <- thisCov < (1-alpha/2) 
  #   
  #   if(lowdef | updef){
  #     
  #     if(lowdef){ 
  #       
  #       lower <- limits[limitX, 1]
  #       upper <- thisQfun(p = 1 - alpha + massLower)
  #     
  #     }else{
  #       
  #       lower <- thisQfun(p = alpha - (1 - massUpper))
  #       upper <- limits[limitX, 2]
  #       
  #     }
  #     
  #     return(Intervals(c(lower, upper)))
  #     
  #   }else{
      
      # grid search adapted from the selectiveInference package
  
  limitsX <- which(limits[,2] > meanEst & limits[,1] < meanEst)
  
  xg = seq(-100*sd, 100*sd, length = gridpts)
  fun = function(x) { 
    tnorm.surv(z = meanEst, mean = x, sd = sd, 
               a = limits[limitsX, 1],
               b = limits[limitsX, 2]) }
  
  int = grid.search(xg, fun, alpha/2, 1-alpha/2, 
                    gridpts, griddepth)
  
  # if(sign(meanEst)==-1) int <- -1*rev(int)
  
  return(Intervals(int))
      
  #   }
  #   
  # }else{
  #   
  #   warning("Can not create one connected confidence interval with ", 1-alpha, " coverage property.\n", 
  #           "Returned interval has ", round(thisCov, 4), " coverage.")
  #   return(limits[limitX,])
  #   
  # }
}


combine_limitObjects <- function(listOfLimitObjects, y){
  
  names <- unique(c(sapply(listOfLimitObjects, names)))
  res <- lapply(names, function(nam){
    
    limits <- lapply(listOfLimitObjects, function(x) x[[nam]])
    
    vT <- limits[[1]]$vT

    limits <- do.call("interval_intersection", lapply(limits, "[[", "limits"))
    
    if(nrow(vT)==1 && !any(sapply(1:nrow(limits),function(i)
      xinInt(x = as.numeric(vT%*%y), int = limits[i,]))))
      stop("Wrong limits. No interval does include the actual value of interest.")
    
    return(list(vT = vT, limits = limits))
    
  })
  names(res) <- names
  
  return(res)
  
}
