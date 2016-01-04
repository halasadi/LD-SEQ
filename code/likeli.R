ll <- function(p, N, n){
  eps = 0.0
  return(sum(log( (1-p)^2*eps^(n) + 2^(1-N)*p*(1-p) + (eps^(N-n))*p*p)))
}


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
mles <- c()
weighted_means <- c()
for (i in 1:2000){
  data <- simulate_data(nind, ptrue, avg_coverage)
  mles[i] <- optimize(f=ll, c(0,1), tol = 1e-5, N=data$N, n = data$n, maximum=TRUE)$maximum
  weighted_means[i] <- sum(data$n)/sum(data$N)
  
}

hist(mles, 20,  main = paste0("emp. mean: ", round(mean(mles), 3)), freq = TRUE)
abline(v = ptrue, col = "red", lwd = 3)
hist(weighted_means, 20, main = paste0("emp. mean: ", round(mean(weighted_means),3)))
abline(v = ptrue, col = "red", lwd = 3)



mse_mle <- mean((mles-ptrue)^2)
print(mse_mle)
mse_wm <-  mean((weighted_means-ptrue)^2)
print(mse_wm)


#### plotting the likleihood surface ####
#p <- seq(0.001, 0.999, 0.01)
#f = rep(0, length(p))
#data <- simulate_data(nind, ptrue, coverage)
#for (i in 1:length(p)){
#  f[i] = ll(p[i], data$N, data$n)
#}
#plot(p,f)
#abline(v = ptrue, col = "red", lwd = 3)