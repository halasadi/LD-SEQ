
## generate test data

simulate_data <- function(nind, ptrue, avg_coverage){
  X <- rbinom(n=nind, size=2, ptrue)
  N <-rpois(nind, avg_coverage)
  n <-rbinom(n=nind, size = N, X/2)  
  # n is the number of non-ref alleles and N is coverage at that site
  ret <- data.frame(N = N, n =n)
  return(ret)
}

ptrue = 0.1
nind <- 5
avg_coverage <- 5
yobs <- c()
epsilon <- c()
nSNPs <- 10
for (i in 1:nSNPs){
  data <- simulate_data(nind, ptrue, avg_coverage)
  yobs[i] <- sum(data$n)/sum(data$N)
  epsilon[i] <- yobs[i]*(1-yobs[i])*sum(data$n+(data$n)^2)/(2*sum(data$n)^2)
}

estimate_sigmabar <- function

