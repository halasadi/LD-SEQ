create_haps <- function(nsamp=100, per){
  # the only type of haplotypes in the pool are either 1-0 or 0-1 (for simplicity)
  # "per" specifies the percent of 1-0 haplotypes  
  haps = matrix(nrow = nsamp, ncol = 2, 0)
  no_10_haps = floor(per*nsamp)
  for (i in 1:no_10_haps){
    haps[i,1] = 1
    haps[i,2] = 0
  }
  for (i in (no_10_haps+1):nsamp){
    haps[i,1] = 0
    haps[i,2] = 1
  }
  return(haps)
}


# potential problem: sigma2 estimated to be < 1 (and even <0 in some cases)
findMLE <- function(obs, mu, Sigma, D){
  
  ll <- function(sigma2){
    cov_m = sigma2*Sigma + D
    sinv = solve(cov_m)
    
    # (from hand-out and verified with Wikipedia)
    logl = log(det(cov_m)) + (t(obs-mu)%*%sinv%*%(obs-mu))
    
    return(as.numeric(logl))
  }
  
  # note: optim actually minimizes
  return(optim(ll, par = 1, method = "L-BFGS-B", lower = 1, upper = 10)$par)
}

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

perform_LDSP <- function(y_obs, nloci, nsamp, haps, pos, n){
  
  # calculate distances between SNPs
  # for example if two SNPs at position 50 and 100
  # then d = 50 
  d = rep(0, nloci-1)
  for (i in 1:(nloci-1)){
    d[i] = as.numeric(pos[i+1])-as.numeric(pos[i])
  }

  epsilon = rep(0, nloci)
  for (i in 1:nloci){
    epsilon[i] = (y_obs[i]*(1-y_obs[i]))/n[i]
  }


  geo_sum = 0
  for (i in 1:(nsamp-1)){
    geo_sum = geo_sum + (1/i)
  }
  geo_sum = (1/geo_sum)
  theta = geo_sum/(nsamp + geo_sum)

  mu = (1-theta)*colMeans(haps) + I(theta/2)

  
  # in this case of only 1_0 and 0_1 haplotypes, cov_panel is singular
  cov_panel = cov(haps)
  
  S = matrix(nrow = nloci, ncol = nloci, 0)

  # rho = 4*N*c_j where c_j is the rate of cross-over per generation per basepair and N is the effective popsize
  # d_j = dist is the physical distance between SNPs j and j+1
  rho = 0 # no recombination
  for (i in 1:nloci){
    for (j in 1:nloci){
      if (i ==j){
        S[i,j] = cov_panel[i,j]
      }
      else{
        # distance between SNPs (e.g. SNP 1 and SNP 10)
        dist = sum(d[min(i,j):(max(i,j)-1)])
        
        # hypothetical effective popsize
        N = 1e4
        
        # since rho = 0, S[i,j] = cov_panel[i,j]
        S[i,j] = exp((-4*N*rho*dist)/nsamp)*cov_panel[i,j]
      }
    }
  }
  
  # no longer singular after the following operation
  Sigma = ((1-theta)^2)*S + I((theta/2)*(1-(theta/2)))
  Sigma[abs(Sigma) < 1e-8] = 0
  
  #sigma2 = findMLE(y_obs, mu, Sigma, diag(epsilon))
  sigma2 = 1
  
  # this is where the dispersion parameter comes in: eqn 13
  d = diag(1/epsilon)
  Sigma_i = solve(Sigma) # there is problems, when nloci is big, kappa(Sigma) is very large (1e36)
  Sigma_bar = solve((1/sigma2)*Sigma_i + d)
  theta_bar = as.vector(Sigma_bar%*%((1/sigma2)*Sigma_i%*%mu + d%*%y_obs))
  
  # posterior
  #return(theta_bar[1])

  mu_t = (mu[1]*Sigma_bar[1,1] - theta_bar[1]*Sigma[1,1])/(Sigma_bar[1,1] - Sigma[1,1])
  
  # I think I missed a minus sign
  sigma_2_t = 1/ ((1/Sigma[1,1]) - 1/(Sigma_bar[1,1]))  
  
  # likelihood + variance
  return(c(mu_t, sigma_2_t))
  
}

mse <- function(y_hat, y){
  return(sum((y_hat-y)^2)/length(y))
}

nloci = 2
# physical position of the two SNPs
pos = c(50, 100)

# starting population
nsamp = 100
# 1_0 haplotype is at 10% frequency and 0_1 is at 90%
haps = create_haps(nsamp, 0.1)

# increase the 1_0 haplotype from 10% to 90% (simulate positive selection)
ev_haps = create_haps(nsamp, 0.9)


## code of testing ##
#haps = matrix(nrow = nsamp, ncol = 2, 0)
#for (i in 1:50){
#  haps[i,1] = 1
#  haps[i,2] = 0
#}
#for (i in 51:60){
#  haps[i,1] = 0
#  haps[i,2] = 1
#}

#for (i in 61:90){
#  haps[i,1] = 1
#  haps[i,2] = 1
#}
#ev_haps = haps
## end testing code ##


lambdas = seq(5, 50, by = 5)
l = length(lambdas)

mse_ldsp_est = rep(0, l)
mse_opt_est = rep(0, l)
mse_obs_est  = rep(0, l)
mean_eff_cov = rep(0, l)

nreps = 1000

for (i in 1:l){
  
  # where to store estimates
  store_ldsp_est = rep(0,nreps)
  store_opt_est = rep(0, nreps)
  store_obs_est = rep(0, nreps)
  store_true_freq = rep(0, nreps)
  store_effective_cov = rep(0, nreps)
  
  for (j in 1:nreps){
    
    # pooling
    lambda = lambdas[i] # coverage
    pooled_info =  pool(lambda, ev_haps, nloci, nsamp)
    n = pooled_info[1,]  # total read counts
    n_1 = pooled_info[2,] # counts of "1" allele
    
    # pseudo-counts
    y_obs = (n_1 + 0.5)/(n + 1)
    
    ldsp_est = perform_LDSP(y_obs, nloci, nsamp, haps, pos, n)
    
    # true pool freqency
    true_freq = colMeans(ev_haps)
    
    # only look at the estimates of SNP 1
    store_true_freq[j] = true_freq[1]
    store_opt_est[j] = (n_1[1] + (n[2]-n_1[2]))/ (n[1] + n[2])
    store_obs_est[j] = y_obs[1]
    #store_obs_est[j] = n_1[1]/n[1]
    store_ldsp_est[j] = ldsp_est[1]
    store_effective_cov[j] = ((1-ldsp_est[1])*ldsp_est[1])/ldsp_est[2]
  }
  
  # calculate MSE
  mse_ldsp_est[i] = mse(store_ldsp_est, store_true_freq)
  mse_opt_est[i] = mse(store_opt_est, store_true_freq)
  mse_obs_est[i] = mse(store_obs_est, store_true_freq)
  mean_eff_cov[i] = mean(store_effective_cov)

}


plot(lambdas, mse_ldsp_est, xlab = "coverage", ylab = "mean square error", col = "red", ylim = c(0,0.05), main = "estimate of SNP 1 frequency", lwd=2.5)
points(lambdas, mse_opt_est, col = "blue", lwd=1.5)
points(lambdas, mse_obs_est, lwd=1.5)
legend("topright", c("read counts only at focal SNP","LDSP", "intuitive optimum"), lty=c(1,1,1), lwd=c(2,2,2),col=c("black","red", "blue"))

# why negative and why not y = 2x?
plot(lambdas, mean_eff_cov, xlab = "coverage", ylab = "effective coverage")