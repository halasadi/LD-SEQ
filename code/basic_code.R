## This script is a pretty much copy and paste of the core componenets of two_loci.R
## but here it's easier to see what's going on

nloci = 2
nsamp = 100

n_10_haps = 70
n_01_haps = nsamp - n_10_haps
haps = c(rep(c(1,0), n_10_haps), c(rep(c(0,1), n_01_haps)))
# put in matrix form
haps = matrix(nrow = nsamp, ncol = 2, haps, byrow=TRUE)

pool <- function(lambda, haps, nloci, nsamp){
  ## simulate pool sequencing ##
  y_obs = rep(0, nloci)
  n = rep(0, nloci)
  n_1 = rep(0, nloci)
  for (i in 1:nloci){
    
    # coverage
    n[i] = rpois(1, lambda)
    
    # ensure number haps sampled is below the number of samples in the pool
    n[i] = min(n[i], nsamp)
    # ensure sample at least 1 haplotype
    n[i] = max(1, n[i])
    
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
theta = geo_sum/(nsamp + geo_sum)
mu = (1-theta)*colMeans(haps) + theta/2

# calculate covariance matrix
cov_panel = matrix(nrow = nloci, ncol = nloci, 0)
f = colMeans(haps)
for (i in 1:nloci){
  for (j in 1:nloci){
    if (i==j){
      cov_panel[i,j] = f[i]*(1-f[i])
    }
    else{
      fij = sum(haps[,i]==1 & haps[,j]==1)/nsamp
      cov_panel[i,j] = fij -f[i]*f[j]
    }
  }
}

S = cov_panel #if rho = 0 then S = cov_panel
Sigma = ((1-theta)^2)*S + (theta/2)*(1-(theta/2))* diag(nloci)

# variances of likelihood
epsilon = y_obs*(1-y_obs)/n

# over-dispersion paramter
sigma2 = 1

# calculate posterior
d = diag(1/epsilon)
Sigma_i = solve(Sigma) 
Sigma_bar = solve((1/sigma2)*Sigma_i + d)
theta_bar = as.vector(Sigma_bar%*%((1/sigma2)*Sigma_i%*%mu + d%*%y_obs))

# likelihood of y_1^true
mu_1 = (mu[1]*Sigma_bar[1,1] - theta_bar[1]*Sigma[1,1])/(Sigma_bar[1,1] - Sigma[1,1])
sigma_2_1 = 1/ ((-1/Sigma[1,1]) + 1/(Sigma_bar[1,1]))  



