#' Forward stepwise AIC selection with adjusted inference for the resulting model
#'
#' @param yname \code{character}; name of response.
#' @param data \code{data.frame}; \code{data.frame} containing response variable and all covariates.
#' @param mod \code{lm} object; smallest linear model.
#' @param var true error variance; if \code{NULL}, the estimate of the final model is used.
#' @param ... further arguments passed to \code{\link{calculate_selinf}}.
#' 
#' @export
#' @importFrom stats update as.formula AIC
#' 
#' @return Returns a list with all visited models and an adjusted inference for the final model 
#' @examples 
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
#' res <- forwardAIC_adjustedInference(yname = "perf",
#'                                     data = cpus,
#'                                     mod = currentmod,
#'                                     var = NULL)
#'                  
#' res$inf
#' 
#' 
#'
forwardAIC_adjustedInference <- function(yname, data, mod, var = NULL, ...)
{
  # get response
  y <- data[[yname]]
  
  # get covariates
  varsInPool <- setdiff(colnames(data), yname)
  
  # start forward stepwise procedure
  improvement <- TRUE
  listOfModelComps <- list() 
  
  while(improvement & length(varsInPool)>0){
    
    # build all other model candidates
    allOtherMods <- lapply(varsInPool, function(thisvar) 
      update(mod, as.formula(paste0(". ~ . + ", thisvar))))
    
    # append mod to list of all model candidates
    listOfModels <- append(allOtherMods, list(mod))
    listOfModelComps <- append(listOfModelComps, list(listOfModels))
    # check AIC
    aics <- sapply(listOfModels, AIC)
    wmaic <- which.min(aics)
    if(wmaic == length(listOfModels)) improvement <- FALSE
    
    # set new best model
    mod <- listOfModels[[wmaic]]
    # get remaining variables
    varsInPool <- varsInPool[-wmaic]
    
  }
  
  # extract test vector
  vTs <- extract_testvec(limo = mod)
  
  # extract components
  compsAIC <- lapply(listOfModelComps, function(lom)
    extract_components(listOfModels = lom,
                       response = y,
                       what = c("aic"))
  )
  
  # calculate limits
  limitsAIC <- lapply(compsAIC, function(comps) calculate_limits(comps, vTs))
  
  sd <- if(is.null(var)) summary(mod)$sigma else sqrt(var)
  
  inf <- calculate_selinf(limitObject = limitsAIC, y = y, sd = sd, ...)
  
  return(list(inf = inf, listOfModels = listOfModelComps))
  
}