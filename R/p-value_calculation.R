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
  
  if((is.na(this_pval) | is.nan(this_pval)) & !lower)
    this_pval <- bryc_tnorm_surv(z = x, sd = sd, 
                                 a = limits[limitX,1], 
                                 b = limits[limitX,2])
  
  if(twosided) this_pval <- 2*min(c(this_pval, 1-this_pval))
  
  return(this_pval)
  
}

#' Calculate confidence interval for truncated normal variable
#' 
#' @param mean mean
#' @param sd standard deviation
#' @param limits \code{Intervals} object containing all limits
#' @param alpha value specifying the interval coverage, which is \code{(1-alpha)}
#' 
#' @examples 
#' library(intervals)
#' # one restriction -- in this case equal to a call to \code{\link{msm::qtnorm}}
#' int1 <- Intervals(c(-3,1))
#' mean = -0.2
#' sd = 1
#' ci_tnorm(mean = mean, sd = sd, limits = int1, alpha = 0.05)
#' 
#' # two intervals -- 1-alpha coverage not possible
#' (int2 <- Intervals(matrix(c(-3,-2,-1,1), 
#'                           ncol = 2, byrow = TRUE)))
#' ci_tnorm(mean = mean, sd = sd, limits = int2, alpha = 0.05)
#' 
#' # two intervals -- 1-alpha coverage possible, but interval not symmetric
#' (int3 <- Intervals(matrix(c(-3,-2.5,-1.9,1), 
#'                           ncol = 2, byrow =  TRUE)))
#' ci_tnorm(mean = mean, sd = sd, limits = int3, alpha = 0.05)
#' 
#' @importFrom msm qtnorm
#' @export
#' 
ci_tnorm <- function(mean, sd, limits, alpha)
{
  
  limits <- limits[order(limits[,1]),]
  pnormstd <- function(x) pnorm((x - mean)/sd)
  limitX <- which(limits[,2] > mean & limits[,1] < mean)
  
  limitXa <- limits[limitX,1]
  massSmallerXa <- if(limitX!=1) 
    sum(sapply(1:(limitX-1), function(i) 
      pnormstd(limits[i,2]) - pnormstd(limits[i,1]))) else 0
  
  nom <- (massSmallerXa + pnormstd(limitXa))
  denom <- sum(sapply(1:nrow(limits), function(i) pnormstd(limits[i,2]) - pnormstd(limits[i,1])))

  massUpper <- (pnormstd(limits[limitX,2]) - nom) / denom
  massLower <- (nom - pnormstd(limits[limitX,1])) / denom
  
  thisQfun <- function(p) qtnorm(p, mean = mean, sd = sd, 
                                 lower = limits[limitX, 1],
                                 upper = limits[limitX, 2])
  thisCov <- massUpper - massLower
  
  ### case discrimination
  
  if(thisCov > 1-alpha){
    
    lowdef <- massLower > (alpha/2)
    updef <- thisCov < (1-alpha/2) 
    
    if(lowdef | updef){
      
      if(lowdef){ 
        
        lower <- limits[limitX, 1]
        upper <- thisQfun(p = 1 - alpha + massLower)
      
      }else{
        
        lower <- thisQfun(p = alpha - (1 - massUpper))
        upper <- limits[limitX, 2]
        
      }
      
      return(Intervals(c(lower, upper)))
      
    }else{
      
      return(Intervals(thisQfun(p = c(alpha/2, 1 - alpha/2))))
      
    }
    
  }else{
    
    warning("Can not create one connected confidence interval with ", 1-alpha, " coverage property.\n", 
            "Returned interval has ", round(thisCov, 4), " coverage.")
    return(limits[limitX,])
    
  }
}

# Adapted from the selectiveInference package
# See https://github.com/selective-inference/R-software/blob/master/selectiveInference
# 
# Returns Prob(Z>z | Z in [a,b]), where mean can be a vector, based on
# A UNIFORM APPROXIMATION TO THE RIGHT NORMAL TAIL INTEGRAL, W Bryc
# Applied Mathematics and Computation
# Volume 127, Issues 23, 15 April 2002, Pages 365--374
# https://math.uc.edu/~brycw/preprint/z-tail/z-tail.pdf

bryc_tnorm_surv <- function(z, mean=0, sd=1, a, b) {
  
  z = (z-mean)/sd
  a = (a-mean)/sd
  b = (b-mean)/sd
  n = length(mean)
  
  ff <- function(z) {
    return((z^2+5.575192695*z+12.7743632)/
             (z^3*sqrt(2*pi)+14.38718147*z*z+31.53531977*z+2*12.77436324))
  }
  
  term1 = exp(z*z)
  o = a > -Inf
  term1[o] = ff(a[o])*exp(-(a[o]^2-z[o]^2)/2)
  term2 = rep(0,n)
  oo = b < Inf
  term2[oo] = ff(b[oo])*exp(-(b[oo]^2-z[oo]^2)/2)
  p = (ff(z)-term2)/(term1-term2)

  p = pmin(1,pmax(0,p))
  return(p)
}



combine_limitObjects <- function(listOfLimitObjects, y){
  
  names <- unique(c(sapply(listOfLimitObjects, names)))
  res <- lapply(names, function(nam){
    
    limits <- lapply(listOfLimitObjects, function(x) x[[nam]])
    
    vT <- limits[[1]]$vT

    limits <- do.call("interval_intersection", lapply(limits, "[[", "limits"))
    
    if(!any(sapply(1:nrow(limits),function(i) 
      xinInt(x = as.numeric(vT%*%y), int = limits[i,])))) 
      stop("Wrong limits. No interval does include the actual value of interest.")
    
    return(list(vT = vT, limits = limits))
    
  })
  names(res) <- names
  
  return(res)
  
}
