ll <- function(ytrue, n, n1){
  eps = 0.0
  return(sum(log( (1-ytrue)^2*eps^(n1) + 2^(1-n)*ytrue*(1-ytrue) + (eps^(n-n1))*ytrue*ytrue)))
}

# MAYBE RE-CHECK FORMULA with MATHEMATICA
fpp <- function(mle, n, n1){
  n0 = n - n1
  a <- sum( (2*(0^n0) + 2*(0^n1) - 2^(2-n)) / (mle^2 * (0^n0) + ((1-mle)^2)*(0^n) + mle*(1-mle)*(2^(1-n))))
  b <- sum( (2*mle*(0^n0) - 2*(1-mle)*(0^n1) + (1-mle)*(2^(1-n)) - mle*(2^(1-n)))^2 / 
              (mle^2 * (0^n0) + ((1-mle)^2)*(0^n1) + mle*(1-mle)*(2^(1-n)))^2)
  return(a-b)
}

LDSEQ <- function(n1, n, Sighat, freq, numSNP){
  # n1    : NOW; vector of length p, FUTURE; dimension n x p, matrix of "1" reads
  # n     : FUTURE; dimension n x p, matrix of total number of reads
  # Sighat: estimated covariance matrix from the panel data using Wen & Stephens, 2011
  # freq  : frequency in the panel
  
  mles <- c()
  epsq <- c()
  tol = 1e-15
  for (j in 1:numSNP){
    mles[j] <- optimize(f=ll, c(0,1), tol = tol, n = n[j], n1 = n1[j], maximum=TRUE)$maximum
    if (mles[j] < tol){
      mles[j] = 0
    }
    epsq[j] <- -1/fpp(mles[j], n = n[j], n1 = n1[j])
  }
  
  ## Crap, it seems the approximation isn't working so well.
  ## the second derivative at the mle is infinite? So var = 1/Inf is 0. 
  
}
