## This script is a pretty much copy and paste of the core componenets of two_loci.R
## but here it's easier to see what's going on
rm(list = ls())
nloci = 2
nsamp = 1000

n_11_haps = 700
n_00_haps = nsamp - n_11_haps
haps = c(rep(c(1,1), n_11_haps), c(rep(c(0,0), n_00_haps)))
# put in matrix form
haps = matrix(nrow = nsamp, ncol = 2, haps, byrow=TRUE)

pool <- function(lambda, haps, nloci, nsamp){
  ## simulate pool sequencing ##
  y_obs = rep(0, nloci)
  n = rep(0, nloci)
  n_1 = rep(0, nloci)
  for (i in 1:nloci){
    # coverage
    n[i] = lambda # rpois(1, lambda)

    # sample locus i from pool
    sampled = sample(haps[,i], n[i])
    
    # number of 1 alleles
    n_1[i] = sum(sampled)
  }
  return(rbind(n, n_1))
}

lambda = 20 # coverage
pooled_info = pool(lambda, haps, 2, nsamp) # perform pooling experiment
n = pooled_info[1,] # total counts
n_1 = pooled_info[2,] # counts of "1" allele
y_obs = (n_1 + 0.5)/(n + 1) # pseudo-counts

geo_sum = 1/sum(1/(1:(nsamp-1)))
#theta = geo_sum/(nsamp + geo_sum)
theta = 0
mu = (1-theta)*colMeans(haps) + theta/2

# calculate panel covariance matrix
cov_panel = cov(haps) * (nsamp-1)/nsamp
S = cov_panel #if rho = 0 then S = cov_panel
Sigma = ((1-theta)^2)*S + (theta/2)*(1-(theta/2))* diag(nloci)

# variances of likelihood
epsilon = y_obs*(1-y_obs)/n

# over-dispersion paramter
sigma2 = 1

# calculate posterior
d = diag(1/epsilon)
Sigma_bar = sigma2*Sigma - sigma2* Sigma %*% solve(diag(1/diag(d)) + sigma2*Sigma) %*% Sigma
#mu_bar2 = solve((1/sigma2)*Sigma_i + d, (1/sigma2)*Sigma_i%*%mu + d%*%y_obs)
mu_bar = solve(diag(nloci) +sigma2 * Sigma %*% d, mu+sigma2*Sigma %*% d %*% y_obs) # system can be singular with real counts

# likelihood y_1^true is normal with following mean and variance
mu_1 = (mu[1]*Sigma_bar[1,1] - mu_bar[1]*Sigma[1,1])/(Sigma_bar[1,1] - Sigma[1,1])
sigma_2_1 = 1/ ((-1/Sigma[1,1]) + 1/(Sigma_bar[1,1]))  

n_e = mu_1*(1-mu_1)/sigma_2_1
print(n_e)
print(n)