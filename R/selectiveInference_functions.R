# Functions copied from selectiveInference package
# See: https://github.com/selective-inference/R-software/
#      blob/master/selectiveInference/R/funs.groupfs.R

TC_surv <- function(TC, sigma, df, E) {
  if (length(E) == 0) {
    stop("Empty TC support")
  }
  
  # Sum truncated cdf over each part of E
  denom <- do.call(sum, lapply(1:nrow(E), function(v) {
    tchi_interval(E[v,1], E[v,2], sigma, df)
  }))
  
  # Sum truncated cdf from observed value to max of
  # truncation region
  numer <- do.call(sum, lapply(1:nrow(E), function(v) {
    lower <- E[v,1]
    upper <- E[v,2]
    if (upper > TC) {
      # Observed value is left of this interval's right endpoint
      if (lower < TC) {
        # Observed value is in this interval
        return(tchi_interval(TC, upper, sigma, df))
      } else {
        # Observed value is not in this interval
        return(tchi_interval(lower, upper, sigma, df))
      }
    } else {
      # Observed value is right of this entire interval
      return(0)
    }
  }))
  
  # Survival function
  value <- numer/denom
  # Force p-value to lie in the [0,1] interval
  # in case of numerical issues
  value <- max(0, min(1, value))
  value
}

tchi_interval <- function(lower, upper, sigma, df) {
  a <- (lower/sigma)^2
  b <- (upper/sigma)^2
  if (b == Inf) {
    integral <- pchisq(a, df, lower.tail = FALSE)
  } else {
    integral <- pchisq(b, df) - pchisq(a, df)
  }
  if ((integral < .Machine$double.eps) && (b < Inf)) {
    integral <- num_int_chi(a, b, df)
  }
  return(integral)
}

num_int_chi <- function(a, b, df, nsamp = 10000) {
  grid <- seq(from=a, to=b, length.out=nsamp)
  integrand <- dchisq(grid, df)
  return((b-a)*mean(integrand))
}

grid.search <- function(grid, fun, val1, val2, gridpts=100, griddepth=2) {
  n = length(grid)
  vals = fun(grid)
  
  ii = which(vals >= val1)
  jj = which(vals <= val2)
  if (length(ii)==0) return(c(grid[n],Inf))   # All vals < val1
  if (length(jj)==0) return(c(-Inf,grid[1]))  # All vals > val2
  # RJT: the above logic is correct ... but for simplicity, instead,
  # we could just return c(-Inf,Inf) 
  
  i1 = min(ii); i2 = max(jj)
  if (i1==1) lo = -Inf
  else lo = grid.bsearch(grid[i1-1],grid[i1],fun,val1,gridpts,
                         griddepth-1,below=TRUE)
  if (i2==n) hi = Inf
  else hi = grid.bsearch(grid[i2],grid[i2+1],fun,val2,gridpts,
                         griddepth-1,below=FALSE)
  return(c(lo,hi))
}

grid.bsearch <- function(left, right, fun, val, gridpts=100, griddepth=1, below=TRUE) {
  n = gridpts
  depth = 1
  
  while (depth <= griddepth) {
    grid = seq(left,right,length=n)
    vals = fun(grid)
    
    if (below) {
      ii = which(vals >= val)
      if (length(ii)==0) return(grid[n])   # All vals < val (shouldn't happen)
      if ((i0=min(ii))==1) return(grid[1]) # All vals > val (shouldn't happen)
      left = grid[i0-1]
      right = grid[i0]
    }
    
    else {
      ii = which(vals <= val)
      if (length(ii)==0) return(grid[1])   # All vals > val (shouldn't happen)
      if ((i0=max(ii))==n) return(grid[n]) # All vals < val (shouldn't happen)
      left = grid[i0]
      right = grid[i0+1]
    }
    
    depth = depth+1
  }
  
  return(ifelse(below, left, right))
}

tnorm.surv <- function(z, mean, sd, a, b, bits=NULL) {
  z = max(min(z,b),a)
  
  # Check silly boundary cases
  p = numeric(length(mean))
  p[mean==-Inf] = 0
  p[mean==Inf] = 1
  
  # Try the multi precision floating point calculation first
  o = is.finite(mean)
  mm = mean[o]
  pp = mpfr.tnorm.surv(z,mm,sd,a,b,bits) 
  
  # If there are any NAs, then settle for an approximation
  oo = is.na(pp)
  if (any(oo)) pp[oo] = bryc.tnorm.surv(z,mm[oo],sd,a,b)
  
  p[o] = pp
  return(p)
}

bryc.tnorm.surv <- function(z, mean=0, sd=1, a, b) {
  z = (z-mean)/sd
  a = (a-mean)/sd
  b = (b-mean)/sd
  n = length(mean)
  
  term1 = exp(z*z)
  o = a > -Inf
  term1[o] = ff(a[o])*exp(-(a[o]^2-z[o]^2)/2)
  term2 = rep(0,n)
  oo = b < Inf
  term2[oo] = ff(b[oo])*exp(-(b[oo]^2-z[oo]^2)/2)
  p = (ff(z)-term2)/(term1-term2)
  
  # Sometimes the approximation can give wacky p-values,
  # outside of [0,1] ..
  #p[p<0 | p>1] = NA
  p = pmin(1,pmax(0,p))
  return(p)
}

ff <- function(z) {
  return((z^2+5.575192695*z+12.7743632)/
           (z^3*sqrt(2*pi)+14.38718147*z*z+31.53531977*z+2*12.77436324))
}

mpfr.tnorm.surv <- function(z, mean=0, sd=1, a, b, bits=NULL) {
  # If bits is not NULL, then we are supposed to be using Rmpf
  # (note that this was fail if Rmpfr is not installed; but
  # by the time this function is being executed, this should
  # have been properly checked at a higher level; and if Rmpfr
  # is not installed, bits would have been previously set to NULL)
  if (!is.null(bits)) {
    z = Rmpfr::mpfr((z-mean)/sd, precBits=bits)
    a = Rmpfr::mpfr((a-mean)/sd, precBits=bits)
    b = Rmpfr::mpfr((b-mean)/sd, precBits=bits)
    return(as.numeric((Rmpfr::pnorm(b)-Rmpfr::pnorm(z))/
                        (Rmpfr::pnorm(b)-Rmpfr::pnorm(a))))
  }
  
  # Else, just use standard floating point calculations
  z = (z-mean)/sd
  a = (a-mean)/sd
  b = (b-mean)/sd
  return((pnorm(b)-pnorm(z))/(pnorm(b)-pnorm(a)))
}
