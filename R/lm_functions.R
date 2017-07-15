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
  ass <- limo$assign
  rlass <- rle(ass)
  rll <- rlass$lengths
  rlv <- rlass$values + 1
  singles <- rlv[rll==1]
  singlesCol <- (1:ncol(X))[sapply(ass+1, function(i) i %in% singles)]
  facs <- rlv[rll>1]
  facsCol <- lapply(facs, function(id) which(ass+1 == id))
  vT <- lapply(singlesCol, function(j) (solve(crossprod(X))%*%t(X))[j,,drop=F])
  Pg <- lapply(facsCol, function(j){
    
    Xj <- X[,j]
    Xminusj <- X[,-j]
    Pminusj <- Xminusj%*%solve(crossprod(Xminusj))%*%t(Xminusj)
    tildeXj <- (diag(nrow(X)) - Pminusj)%*%Xj
    
    tildeXj%*%solve(crossprod(tildeXj))%*%t(tildeXj)
    
  })
  vecs <- c(vT, Pg) 
  vecs <- vecs[order(c(singles,facs))]
  nam <- attr(limo$terms, "term.labels")
  if("(Intercept)"%in%colnames(X)) nam <- c("(Intercept)", nam)
  names(vecs) <- nam
  return(vecs)
  
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
#' @return \code{list} of components for each model. \code{A} and \code{c}: matrices and scalar of 
#' affine inequality; \code{y}: response vector; \code{bestInd} indicator for best model in \code{listOfModels}.
#' 
#' @examples 
#' set.seed(42)
#' y <- rnorm(10)
#' x <- runif(10)
#' mod1 <- lm(y ~ x)
#' 
#' x2 <- runif(10)
#' mod2 <- lm(y ~ x2)
#' 
#' # AIC comparison
#' AIC(mod1,mod2)
#' 
#' lom <- list(mod1, mod2)
#' 
#' # get components
#' comps <- extract_components(lom, response = y)
#' str(comps,1)
#' t(comps$y)%*%comps$A[[1]]%*%comps$y
#' 
#' 
#' @export
#' 
extract_components <- function(listOfModels,
                               response = NULL,
                               what = c("aic", "bic", "llonly", 
                                        "lrt", "Ftest", "sigHunt"),
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
  ps <- ps*REML # if ML was used, this sets all ps to 0
  
  # define the A calculation function
  lb_A <- function(n, p1, p2, pen1, pen2, Px1, Px2){
    
    diag(Px1) <- diag(Px1)-1
    diag(Px2) <- diag(Px2)-1
    
    Px1 <- -1*Px1
    Px2 <- -1*Px2
    
    return(
      (n-p1) * exp( -(p2 - p1 + pen1 - pen2)/n ) * Px2 - (n-p2) * Px1
    )
  }
  
  # prepare components
  if(what == "aic"){ 
    
    critvals <- do.call("AIC", listOfModels)
    bestInd <- which.min(critvals$AIC)
    pens <- 2*critvals$df
    
  }else if(what == "bic"){
    
    critvals <- do.call("BIC", listOfModels)
    bestInd <- which.min(critvals$BIC)
    pens <- log(length(y))*critvals$df
    
  }else if(what == "llonly"){
    
    critvals <- sapply(listOfModels, logLik)
    bestInd <- which.min(critvals)
    
  }else if(what == "lrt"){
    
    if(length(listOfModels) != 2) stop("Likelihood ratio comparisons only allowed for two competing models.")

    p1 <- ncol(model.matrix(listOfModels[[1]]))
    p2 <- ncol(model.matrix(listOfModels[[2]]))
    
    if(p1 > p2) listOfModels <- rev(listOfModels)
        
    critvals <- sapply(listOfModels, logLik)
    pens <- c( pchisq(q = 1 - alpha, df = df), 0)
    bestInd <- as.numeric(-2 * diff(critvals) >= pchisq(q = 1 - alpha, df = df)) + 1
    
    if(bestInd == 2){
      
      # switch components
      pens <- rev(pens)
      ps <- rev(pens)
      modfits <- rev(modfits)
      
    }
     
  }else if(what == "Ftest"){
    
    if(length(listOfModels) != 2) stop("F-test comparisons only allowed for two competing models.")
    
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
    notBestInd <- which((1:2)!=bestInd)
    
    A <- list(modfits[[bestInd]] + kappa * (diag(n) - modfits[[notBestInd]]) - modfits[[notBestInd]])
    
    if((p1 > p2 & bestInd == 1) | (p1 < p2 & bestInd == 2)) A[[1]] <- -1*A[[1]]
    
    
  }else if(what == "sigHunt"){
    
    
    
  }

  if(what != "Ftest"){
  
    # calculate the A's
    A <- lapply((1:length(listOfModels))[-bestInd], function(i)
      lb_A(n = n, p1 = ps[bestInd], p2 = ps[i], pen1 = pens[bestInd],
           pen2 = pens[i], Px1 = modfits[[bestInd]], Px2 = modfits[[i]])
    )
  
  }
  
  # return object components of inequalities
  return(list(A = A, c = 0, y = y, bestInd = bestInd))

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
  
  Pv <- lapply(vTs, function(vt){
    
    if(nrow(vt)!=1){ 
      
      attr(vt, "type") <- "group"
      return(vt) 
      
    }else{ 
      
      r <- t(vt)%*%vt / as.numeric(vt%*%t(vt))
      attr(r, "type") <- "metric"
      return(r)
    }
  
  })
  
  nrMods <- length(comps$A)
  
  res <- lapply(1:length(vTs), function(j){
    
    # for all test vectors (all parameters)
    
    # calculate limits
    vlims <- lapply(1:nrMods, function(k) getBoundsPen(A = comps$A[[k]],
                                                       pv = Pv[[j]], 
                                                       y = comps$y, 
                                                       vt=vTs[[j]])
    )
    
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