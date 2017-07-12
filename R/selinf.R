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


#' Calculate p-values / confidence intervals after likelihood-based or test-based model selection
#' 
#' @param listOfModels either an object of \code{class} \code{limitObject}, the result of the
#' \code{\link{calculate_limits}} function, or a \code{list} of \code{limitObject}s
#' @param y response vector
#' @param sd standard deviation of error used for the p-value calculation (see details)
#' @param alpha value for \code{(1-alpha)}-confidence interval construction. Defaults to \code{0.05}.
#' 
#' @description Wrapper for selinf. See \code{\link{selinf.limitObject}} for more details.
#' 
#'
#' @export
#' @importFrom stats anova logLik model.matrix pchisq pnorm qf resid
#'
selinf_wrapper <- function(listOfModels, response = NULL, 
                           what = c("aic", "bic", "llonly", 
                                    "lrt", "Ftest"), 
                           REML = FALSE, df = 1, 
                           sd, alpha = 0.05)
{
  
  stopifnot(is.list(listOfModels))
  what <- match.arg(what)
  
  comps <- extract_components(listOfModels = listOfModels,
                              response = response,
                              what = what, REML = REML, 
                              alpha = alpha, df =df)
  
  vTs <- extract_testvec(listOfModels[[comps$bestInd]])
  limits <- calculate_limits(comps, vTs)
  
  selinf(limitObject = limits, y = response, sd = sd, alpha = alpha)
  
}
