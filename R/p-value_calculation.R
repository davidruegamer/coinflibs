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
#' @param x quantile
#' @param mean mean
#' @param sd standard deviation
#' @param limits \code{Intervals} object containing all limits
#' @param alpha value specifying the interval coverage, which is \code{(1-alpha)}
#' 
#' @examples 
#' # one restriction
#' int1 <- Intervals(c(-3,1))
#' mean = -0.2
#' sd = 
#' 
#' @importFrom msm qtnorm
ci_tnorm <- function(mean, sd, limits, alpha)
{
  
  limits <- limits[order(limits[,1]),]
  pnormstd <- function(x) pnorm((x - mean)/sd)
  limitX <- which(limits[,2] > mean & limits[,1] < mean)
  
  denom <- sum(sapply(1:nrow(limits), function(i) pnormstd(limits[i,2]) - pnormstd(limits[i,1])))
  massUpper <- pnormstd(limits[limitX,2])
  massLower <- pnormstd(limits[limitX,1])
  
  thisQfun <- function(p) qtnorm(p, mean = mean, sd = sd, 
                                 lower = limits[limitX, 1],
                                 upper = limits[limitX, 2])
  thisCov <- (massUpper - massLower) / denom
  
  ### case discrimination
  
  if(thisCov > 1-alpha){
    
    lowdef <- massLower > (alpha/2)
    updef <- massUpper < (1-alpha/2) 
    
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
            "Returned interval has ", thisCov, " coverage.")
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

#' Calculate p-values / confidence intervals after likelihood-based or test-based model selection
#' 
#' @param limitObject either an object of \code{class} \code{limitObject}, the result of the
#' \code{\link{calculate_limits}} function, or a \code{list} of \code{limitObject}s
#' @param y response vector
#' @param sd standard deviation of error used for the p-value calculation (see details)
#' @param alpha value for \code{(1-alpha)}-confidence interval construction. Defaults to \code{0.05}.
#' 
#' @description This function takes an \code{limitObject}, which is produced by the 
#' \code{\link{calculate_limits}} (or a list of \code{limitObject}s) and calculates
#' the selective p-value / confidence intervals based on the given limitation and the true 
#' residual standard deviation \code{sd}. Since the true standard deviation is usually unknown in practice
#' one can plug in an estimate for the true standard deviation. This approach, however,
#' strongly depends on the goodness of the estimate and in particular produces unreliable results 
#' if the number of covariates is relatively large in comparison to the number of observations.
#' 
#' @examples 
#'
#' 
#' library(MASS)
#' # Fit initial model
#' cpus$perf <- log10(cpus$perf)
#' mod <- lm(perf ~ .-name, data = cpus)
#' 
#' # use the stepAIC function to find the best model in a backward 
#' # stepwise search
#' cpus.lm <- stepAIC(mod, trace = FALSE, direction = "backward")
#' 
#' # recalculate all visited models
#' lom <- c(lapply(colnames(model.matrix(mod)), function(x) 
#'   update(mod, as.formula(paste0("perf ~ .-", x)))), list(mod))
#'   
#' # extract the components of all visited models
#' compsAIC <- extract_components(listOfModels = lom,
#'                                response = cpus$perf,
#'                                what = c("aic"))
#'                                  
#' # perform likelihood ratio test at level 
#' alpha = 0.001
#' 
#' # check for non-significant variables
#' coefTable <- summary(cpus.lm)$coefficients
#' drop <- rownames(coefTable)[alpha < coefTable[,4]]
#' 
#' # drop non-significant variable
#' cpus.lm2 <- update(cpus.lm, as.formula(paste0(".~.-",drop)))
#' 
#' # extract components associated with the LRT comparison
#' compsLRT <- extract_components(list(cpus.lm, cpus.lm2),
#'                                response = cpus$perf, 
#'                                what = "lrt",
#'                                alpha = alpha)
#'                                
#' # naive inference
#' unadj_pvs <- summary(cpus.lm2)$coefficients[,4]
#' 
#' # now extract testvector, calculate limits and perform selective tests
#' # test vectors
#' 
#' vTs <- extract_testvec(cpus.lm2)
#' 
#' # calculate limits
#' limitsAIC <- calculate_limits(compsAIC, vTs)
#' limitsLRT <- calculate_limits(compsLRT, vTs)
#' 
#' # check restriction on p-values separately
#' cbind(
#' selinf(limitObject = limitsAIC, y = cpus$perf, sd = sigma(cpus.lm2)),
#' unadjusted_pval = unadj_pvs
#' )
#' 
#' cbind(
#' selinf(limitObject = limitsLRT, y = cpus$perf, sd = sigma(cpus.lm2)),
#' unadjusted_pval = unadj_pvs
#' )
#' 
#' # calculate p-values (does automatically combine limitObjects)
#' res <- selinf(limitObject = list(limitsAIC, limitsLRT), 
#'                        y = cpus$perf, 
#'                        # plugin estimate for true value
#'                        sd = sigma(cpus.lm2))
#'                        
#' cbind(res, unadjusted_pval = unadj_pvs)
#'
#'
#' @export
#' @importFrom stats anova logLik model.matrix pchisq pnorm qf resid
#'
selinf <- function(limitObject, y, sd, alpha = 0.05)
{
  
  if(class(limitObject) != "limitObject"){
    
    
    ### list of limitObjects -> combine limits
    ### stop if that's not the case
    if(any(sapply(limitObject, class) != "limitObject"))
      stop("limitObject must either be of class 'limitObject' or a list with 'limitObject's")
    
    limitObject <- combine_limitObjects(limitObject, y = y)
    
  }
  
  pvals <- lowLim <- upLim <- ests <- rep(NA, length(limitObject))
  
  # calculate p-values for all limitObjects
  for(j in 1:length(limitObject)){
    
    vT <- limitObject[[j]]$vT
    mu <- as.numeric(vT%*%y)
    this_sd <- sd * sqrt(as.numeric(tcrossprod(vT)))
    
    limits <- limitObject[[j]]$limits
    
    pvals[j] <- mult_tnorm_surv(x = mu, mean = 0, 
                                sd = this_sd,
                                limits = limits)
    
    int <- ci_tnorm(mean = mu, sd = this_sd, 
                    limits = limits, alpha = alpha)
    
    lowLim[j] <- int[[1]]
    upLim[j] <- int[[2]]
    
    ests[j] <- mu
    
  }
  
  return(
    data.frame(varname = names(limitObject),
               teststat = ests,
               lower = lowLim,
               upper = upLim,
               pval = pvals
    ))
  
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
