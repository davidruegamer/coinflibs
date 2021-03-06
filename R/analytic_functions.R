xinInt <- function(x,int) any(int[,2] > x & int[,1] < x)

solveQuadIneq <- function(A, c, Pv, PvO = NULL, y)
{
  
  if(is.null(PvO)) PvO <- diag(length(y)) - Pv
  alpha <- as.numeric(t(y)%*%Pv%*%A%*%Pv%*%y)
  beta <- as.numeric(2*t(y)%*%Pv%*%A%*%PvO%*%y)
  gamma <- as.numeric(t(y)%*%PvO%*%A%*%PvO%*%y + c)
  
  dis <- beta^2 - 4*alpha*gamma
  if(dis < 0) return(c(-Inf,Inf))
  
  sort(c(
    (-beta - sqrt(dis))/(2*alpha),
    (-beta + sqrt(dis))/(2*alpha)
  ))
  
}

# solveQuadIneqGen <- function(A, c, gammaj, y)
# {
#   
#   gAg <- as.numeric(t(gammaj)%*%A%*%gammaj)
#   yAy <- as.numeric(t(y)%*%A%*%y)
#   
#   alpha <- -1*gAg
#   beta <- 2*gAg
#   gamma <- yAy - gAg + c
#   
#   det <- beta^2 - 4*alpha*gamma
#   if(det < 0) return(c(-Inf,Inf))
#   
#   sort(c(
#     (-beta - sqrt(det))/(2*alpha),
#     (-beta + sqrt(det))/(2*alpha)
#   ))
#   
# }




# getBoundsPen <- function(sigmasq1, sigmasq2, 
#                          Px1, Px2, pen1, pen2, 
#                          pv, y, vt, REML)
# {
#   
#   
#   # a <- (sigmasq1-sigmasq2)
#   # 
#   # b <- 2*(sigmasq2*(mu1) - sigmasq1*(mu2))
#   # 
#   # c <- ((n*log(sigmasq2/sigmasq1) + pen2 - pen1)*(sigmasq1*sigmasq2) -
#   #         sigmasq2*as.numeric(crossprod(mu1)) + sigmasq1*as.numeric(crossprod(mu2)))
#   
#   n <- length(y)
#   
#   A <- sigmasq1 * (diag(n) - Px2) - sigmasq2 * (diag(n) - Px1)
#   
#   c <- (n*log(sigmasq2/sigmasq1) + pen2 - pen1)*(sigmasq1*sigmasq2)
#   
#   # stopifnot(as.numeric(t(y)%*%A%*%y + c)>=0)
#   
#   taus <- solveQuadIneq(A=A, c=c, Pv=pv, y=y)
#   
#   mu <- as.numeric(vt%*%y)
#   
#   Vlo <- mu*taus[1]
#   Vup <- mu*taus[2]
#   
#   stopifnot(length(taus)>0)
#   
#   if(mu < 0){
#     
#     tt <- Vlo
#     Vlo <- Vup
#     Vup <- tt
#     
#   }
#   
#   if((mu > Vup & Vup > Vlo) | (mu < Vlo & Vlo < Vup)){
#     
#     resVec <- Intervals(rbind(c(-Inf, Vlo), c(Vup, Inf)))
#     type <- "dyadic"
#     
#   }else if(Vlo < mu & mu < Vup){
#     
#     resVec <- Intervals(c(Vlo,Vup))
#     type <- "single"
#     
#   }else if(Vlo > Vup){
#     
#     stop("Something went wrong.")
#     
#   }
#   
#   if(!xinInt(x = mu, int = resVec)) stop("Something went wrong.")
#   
#   return(list(lim=resVec,type=type))
#   
#   
# }

#' Function to calculation boundaries for given components
#' 
#' @param A restriction matrix in affine inequality
#' @param pv projection matrix for test vector or its left singular values for group variables
#' @param y response
#' @param vt test vector
#' @param c additive value in affine inequality; defaults to 0
#' 
#' @importFrom intervals Intervals
#' 
#'
getBoundsPen <- function(A, pv, y, vt, c = 0)
{
  
  n <- length(y)
  
  if(attr(pv, "type") == "group"){
    
    Ugtilde <- svd(pv)$u
    pv <- tcrossprod(Ugtilde)
    pvo <- diag(n) - pv
    R <- t(Ugtilde) %*% y
    TC <- sqrt(sum(R^2))
    pv <- pv / TC
    
    
  }else{
    
    pvo <- NULL
    
  }
  
  taus <- solveQuadIneq(A = A, c = c, Pv = pv, PvO = pvo, y = y)

  if(!is.null(pvo)){ 
    
    mu <- TC
    Vlo <- taus[1]
    Vup <- taus[2]
    
  }else{
    
    mu <- as.numeric(vt%*%y)
    Vlo <- mu*taus[1]
    Vup <- mu*taus[2]
    
  }
  
  stopifnot(length(taus)>0)
  
  if(mu < 0){
    
    tt <- Vlo
    Vlo <- Vup
    Vup <- tt
    
  }
  
  if((mu > Vup & Vup > Vlo) | (mu < Vlo & Vlo < Vup)){
    
    resVec <- Intervals(rbind(c(-Inf, Vlo), c(Vup, Inf)))
    type <- "dyadic"
    
  }else if(Vlo < mu & mu < Vup){
    
    resVec <- Intervals(c(Vlo,Vup))
    type <- "single"
    
  }else if(Vlo > Vup){
    
    stop("Something went wrong.")
    
  }
  
  if(!xinInt(x = mu, int = resVec)) stop("Something went wrong.")
  
  return(list(lim = resVec, type = type))
  
  
}
