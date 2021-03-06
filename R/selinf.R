#' Calculate p-values / confidence intervals after likelihood-based or test-based model selection
#' 
#' @param limitObject either an object of \code{class} \code{limitObject} (the result of the
#' \code{\link{calculate_limits}} function) or a \code{list} of \code{limitObject}s
#' @param y response vector
#' @param sd standard deviation of error used for the p-value calculation (see details)
#' @param alpha value for \code{(1-alpha)}-confidence interval construction. Defaults to \code{0.05}.
#' @param ... options passed to \code{\link{ci_tnorm}}.
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
#' calculate_selinf(limitObject = limitsAIC, y = cpus$perf, 
#' sd = summary(cpus.lm2)$sigma),
#' unadjusted_pval = unadj_pvs
#' )
#' 
#' cbind(
#' calculate_selinf(limitObject = limitsLRT, y = cpus$perf, 
#' sd = summary(cpus.lm2)$sigma),
#' unadjusted_pval = unadj_pvs
#' )
#' 
#' # calculate p-values (does automatically combine limitObjects)
#' res <- calculate_selinf(limitObject = list(limitsAIC, limitsLRT), 
#'                        y = cpus$perf, 
#'                        # plugin estimate for true value
#'                        sd = summary(cpus.lm2)$sigma)
#'                        
#' cbind(res, unadjusted_pval = unadj_pvs)
#'
#'
#'#########################################################
#'### Another example for a forward stepwise AIC search ###
#'#########################################################
#'
#' library(MASS)
#' # use the cpus data
#' data("cpus")
#'
#' # Fit initial model
#' cpus$perf <- log10(cpus$perf)
#' cpus$cach <- as.factor(cpus$cach)
#' cpus$name <- NULL
#' currentmod <- lm(perf ~ 1, data = cpus)
#'
#' # make a stepwise AIC-based forward search
#' # for all variables in the pool of possible covariates
#' varsInPool <- colnames(cpus)[-7]
#'
#' # since the stepAIC function does not provide the models
#' # fitted in each step, we have to do the search 'manually'
#' improvement <- TRUE
#' listOfModelComps <- list() 
#'
#' # do the forward stepwise AIC search...
#' while(improvement & length(varsInPool)>0){
#'
#'  # compute all other models
#'  allOtherMods <- lapply(varsInPool, function(thisvar) 
#'                         update(currentmod, as.formula(paste0(". ~ . + ", thisvar))))
#'                         
#'  # store all models that were examined in this step  
#'  listOfModels <- append(allOtherMods, list(currentmod))
#'  # save this list for later
#'  listOfModelComps <- append(listOfModelComps, list(listOfModels))
#'  
#'  # check the AIC of all models
#'  aics <- sapply(listOfModels, AIC)
#'  # what is the best model?
#'  (wmaic <- which.min(aics))
#'  # is there any improvement?
#'  if(wmaic == length(listOfModels)) improvement <- FALSE
#'  # redefine the current (best) model
#'  currentmod <- listOfModels[[wmaic]]
#'  # and update the variables available
#'  varsInPool <- varsInPool[-wmaic]
#'
#' }
#'
#' # variables left, which did not improve the model
#' varsInPool
#' # the final model call
#' currentmod$call
#'
#' # get the test vector from the current model
#' vTs <- extract_testvec(limo = currentmod)
#'
#' # extract list of model components in each step when comparisons
#' # are done based on the AIC
#' listOfComps <- lapply(listOfModelComps, function(lom) 
#'  extract_components(listOfModels = lom, response = cpus$perf, what = "aic"))
#'  
#' # calculate the truncation limits for each of the comparisons in each iteration
#' listOfLimits <- lapply(listOfComps, function(lom) 
#'   calculate_limits(comps = lom, vTs = vTs))
#'   
#' # now compute selective inference, which adjust for the forward stepwise AIC search 
#' # by supplying the lists of limits
#' calculate_selinf(limitObject = listOfLimits,
#'                  y = cpus$perf, 
#'                  sd = sigma(currentmod))
#'                  
#' ########################################################
#' # now do that with the function provided in the package
#' ########################################################
#' 
#' #' library(MASS)
#' # use the cpus data
#' data("cpus")
#'
#' # Fit initial model
#' cpus$perf <- log10(cpus$perf)
#' cpus$cach <- as.factor(cpus$cach)
#' cpus$name <- NULL
#' currentmod <- lm(perf ~ 1, data = cpus)
#' 
#' res <- forwardAIC_adjustedInference(yname = "perf",
#'                                     data = cpus,
#'                                     mod = currentmod,
#'                                     var = NULL)
#'                  
#' res$inf
#'
#' @export
#' @importFrom stats anova logLik model.matrix pchisq pnorm qf resid
#'
calculate_selinf <- function(limitObject, y, sd, alpha = 0.05, ...)
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
    limits <- limitObject[[j]]$limits
    
    
    if(nrow(vT)==1){ # truncated normal case
      
      mu <- as.numeric(vT%*%y)
      this_sd <- sd * sqrt(as.numeric(tcrossprod(vT)))
      
      pvals[j] <- mult_tnorm_surv(x = mu, mean = 0, 
                                  sd = this_sd,
                                  limits = limits)
      
      int <- ci_tnorm(meanEst = mu, sd = this_sd, 
                      limits = limits, alpha = alpha, ...)
      
      lowLim[j] <- int[[1]]
      upLim[j] <- int[[2]]
      
      ests[j] <- mu
    
    }else{
  
      Ugtilde <- svd(vT)$u
      R <- t(Ugtilde) %*% y
      
      lowLim[j] <- upLim[j] <- NA # tbd
      ests[j] <- sqrt(sum(R^2))
      
      pvals[j] <- TC_surv(TC = ests[j], sigma = sd, 
                          df = ncol(Ugtilde), E = limits)
      
      
    } 
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
#' @param listOfModels a \code{list} of models representing one comparison
#' @param ... optionally more \code{list}s of models representing further comparisons
#' @param response response vector
#' @param what \code{character} describing the method, on the basis of which variable selection was performed.
#' Supported variable selection procedures are AIC (\code{"aic"}), BIC (\code{"bic"}), 
#' likelihood-based selection (\code{"llonly"}) and significance tests based on the likelihood-ratio 
#' or the F-statistic (\code{"lrt"}, \code{"Ftest"}). If different selection procedures have been used, 
#' please supply the corresponding character for each list of models.
#' @param REML \code{logical}; whether scale estimates should be calculated via \code{REML} or not. Default is
#' \code{FALSE} since the log-likelihood is calculated with \code{REML=FALSE} by the \code{logLik} function.
#' @param df positive integer; defines the degree of freedom when using the likelihood ratio test for model selection.
#' Defaults to \code{1} (testing one parameter at a time).
#' @param sd standard deviation of error used for the p-value calculation (see details)
#' @param alpha value for \code{(1-alpha)}-confidence interval construction. Defaults to \code{0.05}.
#' @param gridpts,griddepth,range options for the calculation of confidence intervals. 
#' See \code{\link{ci_tnorm}} for more details.
#' 
#' @description This is a wrapper and convenience function for several functions included in this package. 
#' Please see \code{\link{calculate_selinf}} for more details and how to calculate inference when
#' several model selection procedures are combined. 
#' When more than one list is supplied, the best model in the first \code{listOfModels} is used to perform inference.
#' 
#' 
#'
#' @export
#' @importFrom stats anova logLik model.matrix pchisq pnorm qf resid dchisq
#' @examples 
#' library(MASS)
#' data("cpus")
#' # Fit initial model
#' cpus$perf <- log10(cpus$perf)
#' cpus$cach <- as.factor(cpus$cach)
#' mod <- lm(perf ~ .-name, data = cpus)
#' 
#' # use the stepAIC function to find the best model in a backward 
#' # stepwise search
#' cpus.lm <- stepAIC(mod, trace = FALSE, direction = "backward", steps = 3)
#' # check model selection
#' cpus.lm$anova$Step
#' 
#' # recalculate all visited models in the first step
#' lom1 <- c(lapply(attr(mod$terms, "term.labels"), function(x) 
#' update(mod, as.formula(paste0("perf ~ .-", x)))), list(mod))
#'
#' # perform likelihood ratio test at level 
#' alpha = 0.001
#' 
#' # check for non-significant variables
#' coefTable <- anova(cpus.lm)
#' drop <- rownames(coefTable)[alpha < coefTable[-nrow(coefTable),5]]
#' 
#' # drop non-significant variable
#' cpus.lm2 <- update(cpus.lm, as.formula(paste0(".~.-",drop)))
#'
#' # compute selective inference
#' selinf(list(cpus.lm, cpus.lm2), lom1, 
#' response = cpus$perf,
#' what = c("Ftest", "aic"),
#' sd = summary(cpus.lm2)$sigma)
#'
selinf <- function(listOfModels, ...,
                   response = NULL, 
                   what,
                   REML = FALSE, df = 1, 
                   sd, alpha = 0.05, 
                   gridpts = 500, griddepth = 3, range = 100)
{
  
  stopifnot(is.list(listOfModels))
  
  if(!missing(...))
  {
    
    comps <- extract_components(listOfModels = listOfModels,
                                response = response,
                                what = what[[1]], REML = REML, 
                                alpha = alpha, df = df)  
    
    vTs <- extract_testvec(listOfModels[[comps$bestInd]])
    
    lim1 <- calculate_limits(comps, vTs)

    ll <- list(...)
    if(length(what)>1){
      if(length(what)!=length(ll)+1)
      {
        
        stop(paste0("If 'what' is different for each list of models, ",
                    "it must be a vector of the same length as supplied lists.")
        )
        
      }else{
        
        # selection procedure is the same for each comparison
        # -> replicate what
        what <- rep(what, length(ll))
      
      }
    }
    
    limits <- lapply(1:length(ll), function(i){
      
      listItem <- ll[[i]]
      
      comps <- extract_components(listOfModels = listItem,
                                  response = response,
                                  what = what[[i+1]], 
                                  REML = REML, 
                                  alpha = alpha, df = df)
      
      calculate_limits(comps, vTs)
      
    })
    
    limits <- c(list(lim1), limits)
    
  }else{
  
    comps <- extract_components(listOfModels = listOfModels,
                                response = response,
                                what = what, REML = REML, 
                                alpha = alpha, df = df)
    
    vTs <- extract_testvec(listOfModels[[comps$bestInd]])
    limits <- calculate_limits(comps, vTs)
    
  }
  
  calculate_selinf(limitObject = limits, y = response, 
                   sd = sd, alpha = alpha,
                   gridpts = gridpts, griddepth = griddepth, 
                   range = range)
  
}
