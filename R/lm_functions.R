#' extract all test vectors from \code{lm} object
#' 
#' @param limo \code{lm} object
#' @return \code{list} of test vectors
#' @export
#' 
#' @examples
#' set.seed(42)
#' y <- rnorm(10)
#' x <- runif(10)
#' mod <- lm(y ~ x)
#' 
#' # get test vectors
#' vT <- extract_testvec(mod)
#' # calculate coefficient estimates
#' sapply(vT, function(vt) vt%*%y)
#' coef(mod)
#' 
extract_testvec <- function(limo)
{
  
  X <- model.matrix(limo)
  vT <- lapply(1:(ncol(X)), function(j) (solve(crossprod(X))%*%t(X))[j,,drop=F])
  names(vT) <- colnames(X)
  return(vT)
  
}

#' hat matrix calculation of linear model
#' 
#' @param limo \code{lm} object
#' @param X \code{matrix} respresenting the design matrix of a linear model
#' @return hat matrix
#' 
#' @examples
#' set.seed(42)
#' y <- rnorm(10)
#' x <- runif(10)
#' mod <- lm(y ~ x)
#' 
#' str(getHatMat(mod),1)
#' str(getHatMat(X = as.matrix(cbind(interc=rep(1,10),x))),1)
#' 
#' @export
getHatMat <- function(limo = NULL, X = NULL)
{
  
  if(is.null(X)){
    
    if(is.null(limo)) stop("Either limo or X argument must be specified.")
    X <- model.matrix(limo)
    
  } 
  
  if(ncol(X)==0){
    
    # null model
    return(matrix(0, nrow=nrow(X), ncol=nrow(X)))
    
  }else{
    
    tX <- t(X)
    return(X%*%solve(tX%*%X)%*%tX)
    
  }
  
  

}

#' extract components from list of \code{lm} objects
#' 
#' @param listOfModels \code{list} of \code{lm} objects
#' @param response response vector; must be supplied if not stored in the \code{lm} object.
#' @param what \code{character} describing the method, on the basis of which variable selection was performed
#' @param REML \code{logical}; whether scale estimates should be calculated via \code{REML} or not. Default is
#' \code{FALSE} since the log-likelihood is calculated with \code{REML=FALSE} by the \code{logLik} function.
#' @param alpha significance level for selection via likelihood ratio test or significance hunting.
#' @param df positive integer; defines the degree of freedom when using the likelihood ratio test for model selection.
#' Defaults to \code{1} (testing one parameter at a time).
#' @return \code{list} of components for each model
#' @examples 
#' set.seed(42)
#' y <- rnorm(10)
#' x <- runif(10)
#' mod1 <- lm(y ~ x)
#' 
#' x2 <- runif(10)
#' mod2 <- lm(y ~ x2)
#' 
#' AIC(mod1,mod2)
#' 
#' lom <- list(mod1, mod2)
#' 
#' # get components
#' comps <- extract_components(lom, response = y)
#' str(comps,1)
#' 
#' @export
#' 
extract_components <- function(listOfModels,
                               response = NULL,
                               what = c("aic", "bic", "llonly", 
                                        "lrt", "Ftest"),
                               REML = FALSE, alpha = NULL, 
                               df = 1)
{
  
  
  # get response from model
  y <- listOfModels[[1]]$y
  if(is.null(y)){
    if(is.null(response)) stop("response must be supplied") else
      y <- response  
  } 
  # get number of observations
  n <- length(y)
  
  # check alpha and df
  if(!is.null(alpha)) stopifnot(alpha > 0 & alpha <= 0.5)
  if(!is.null(df) && 
     (all.equal(round(df),df)!=TRUE | df < 0)) stop("df must be a positive integer.") 
  
  # check what model selection procedure was conducted
  what <- match.arg(what)
  
  #####################################
  # define components
  
  # get hat matrix of each model
  modfits <- lapply(listOfModels, getHatMat)
  
  # define number of columns per model
  ps <- sapply(listOfModels, function(x) ncol(model.matrix(x)))
  
  if(what == "aic"){ 
    
    critvals <- do.call("AIC", listOfModels)
    bestInd <- which.min(critvals$AIC)
    pens <- 2*critvals$df
    
  }else if(what == "bic"){
    
    critvals <- do.call("BIC", listOfModels)
    bestInd <- which.min(critvals$BIC)
    pens <- log(length(y))*critvals$df
    
  }else if(what == "none"){
    
    critvals <- sapply(listOfModels, logLik)
    bestInd <- which.min(critvals)
    pens <- rep(0, length(critvals))
    
  }else if(what == "lrt"){
    
    stopifnot(length(listOfModels)==2)
    critvals <- sapply(listOfModels, logLik)
    
    pens <- c( pchisq(q = 1 - alpha, df = df), 0)
    
    # if models are given in the order H_0, H_1
    # take the negative LR
    bestInd <- as.numeric(-2 * diff(critvals) >= pchisq(q = 1 - alpha, df = df)) + 1
    
    
  }else if(what == "Ftest"){
    
    stopifnot(length(listOfModels)==2)
    
    p1 <- ncol(model.matrix(listOfModels[[1]]))
    p2 <- ncol(model.matrix(listOfModels[[2]]))
    
    testres <- anova(listOfModels[[1]], listOfModels[[2]])
    critval <- qf(0.95, 1, n-p2)
    kappa <- critval * (p2-p1)/(n-p2)
    # check
    bestInd <- sum(resid(listOfModels[[1]])^2) -
      (1+kappa)*sum(resid(listOfModels[[2]])^2) <= 0
    stopifnot((testres$p.value <= 0.05) == !bestInd)
    
    # if models are given in the order H_0, H_1
    # take the negative LR
    bestInd <- 2 - bestInd
    
    # trick to handle Ftest with the same functions later
    pens <- c(0,-n*log(1/kappa))
    modfits[[1]] <- diag(n) + modfits[[2]] - modfits[[1]]
    
  }

  # return objects
  return(list(ps = ps,
              bestInd = bestInd,
              modfits = modfits,
              pens = pens,
              y = y,
              REML = REML))

}




#' Calculate the space restriction limits of a parameter
#'
#' @param comps result of \code{calculate_components} function
#' @param vTs named list of testvectors
#' @return named \code{list} of class \code{limitObject}. Each element corresponds to one 
#' test vector and does contain \code{vT}, the test vector, and \code{limits}, 
#' an \code{Intervals} object for the calculated limits
#'
#' @importFrom intervals interval_intersection
#' @export
#'
calculate_limits <- function(comps, vTs)
{
  
 
  if(!is.list(vTs)) vTs <- list(vTs)
  
  Pv <- lapply(vTs, function(vt) t(vt)%*%vt / as.numeric(vt%*%t(vt)))
  
  nrMods <- length(comps$ps)
  bestInd <- comps$bestInd
  notbestInd <- (1:nrMods)[-comps$bestInd]
  
  res <- lapply(1:length(vTs), function(j){
    
    # for all test vectors (all parameters)
    
    # calculate limits
    vlims <- lapply(notbestInd, 
               function(k) 
                 getBoundsPen(
                 p1 = comps$ps[bestInd],
                 p2 = comps$ps[k],
                 Px1 = comps$modfits[[bestInd]],
                 Px2 = comps$modfits[[k]],
                 pen1 = comps$pens[bestInd],
                 pen2 = comps$pens[k],
                 pv = Pv[[j]], 
                 y = comps$y, 
                 vt=vTs[[j]], 
                 REML = comps$REML))#)
    
    # intersect limits
    limits <- do.call("interval_intersection", lapply(vlims, "[[", "lim"))
    
    return(
      list(vT = vTs[[j]], 
           limits = limits)
    )
    
  })
  
  names(res) <- names(vTs)
  class(res) <- "limitObject"
  return(res)
  
  
}