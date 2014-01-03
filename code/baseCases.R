create_haps <- function(nsamp=100, per){
  # the only type of haplotypes in the pool are either 1-0 or 0-1 (for simplicity)
  # per specifies the percent of 1-0 haplotypes  
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

# global variables are only to debug code
#global_obs = 0
#global_mu = 0
#global_Sigma = 0
#lobal_D = 0

# PROBLEM HERE: sigma2 estimated to be < 1
findMLE <- function(obs, mu, Sigma, D){
  #global_obs <<- obs
  #global_mu <<- mu
  #global_Sigma <<- Sigma
  #global_D <<- D
  
  ll <- function(sigma2){
    cov_matrix = sigma2*Sigma + D
    sinv = solve(cov_matrix)
    # equation from handout (agrees with Wikipedia)
    logl = log(det(cov_matrix)) + (t(obs-mu)%*%sinv%*%(obs-mu))
    return(as.numeric(logl))
  }
  
  # note: optim actually minimizes
  return(optim(ll, par = 1, method = "L-BFGS-B", lower = 0.1, upper = 10)$par)
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
        
        S[i,j] = exp((-4*N*rho*dist)/nsamp)*cov_panel[i,j]
      }
    }
  }

  Sigma = ((1-theta)^2)*S + I((theta/2)*(1-(theta/2)))
  Sigma[abs(Sigma) < 1e-8] = 0
  
  #sigma2 = findMLE(y_obs, mu, Sigma, diag(epsilon))
  #print(sigma2)
  sigma2 = 1
  
  # this is where the dispersion parameter comes in: eqn 13
  d = diag(1/epsilon)
  Sigma_i = solve(Sigma) # there is problems, when nloci is big, kappa(Sigma) is very large (1e36)
  Sigma_bar = solve((1/sigma2)*Sigma_i + d)
  theta_bar = as.vector(Sigma_bar%*%((1/sigma2)*Sigma_i%*%mu + d%*%y_obs))
  return(as.vector(theta_bar))
}

mse <- function(y_hat, y){
  return(sum((y_hat-y)^2)/length(y))
}

nloci = 2
# physical position of the two SNPs
pos = c(50, 100)

# starting population
nsamp = 100
haps = create_haps(nsamp, 0.4)

# increase the 1-0 haplotype from 40% to 90% (simulate positive selection)
ev_haps = create_haps(nsamp, 0.9)

lambdas = seq(5, 50, by = 5)
l = length(lambdas)

mse_est = rep(0, l)
mse_opt = rep(0, l)
mse_obs  = rep(0, l)

nreps = 1000

for (i in 1:l){
  store_ldsp_est = rep(0,nreps)
  store_opt_est = rep(0, nreps)
  store_obs_est = rep(0, nreps)
  store_true_freq = rep(0, nreps)
  
  for (j in 1:nreps){
    
    # pooling
    lambda = lambdas[i] # coverage
    pooled_info =  pool(lambda, ev_haps, nloci, nsamp)
    n = pooled_info[1,]  # total read counts
    n_1 = pooled_info[2,] # counts of "1" allele
    
    # psedu-counts
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
  }
  
  # calculate MSE
  mse_est[i] = mse(store_ldsp_est, store_true_freq)
  mse_opt[i] = mse(store_opt_est, store_true_freq)
  mse_obs[i] = mse(store_obs_est, store_true_freq)

}


plot(lambdas, mse_est, xlab = "coverage", ylab = "mean square error", col = "red", ylim = c(0,0.03), main = "estimate of SNP 1 frequency", lwd=1.5)
points(lambdas, mse_opt, col = "blue", lwd=1.5)
points(lambdas, mse_obs, lwd=1.5)
legend("topright", c("read counts only at focal SNP","LDSP", "intuitive optimum"), lty=c(1,1,1), lwd=c(2,2,2),col=c("black","red", "blue"))
#legend("topright", c("read counts only at focal SNP","LDSP"), lty=c(1,1), lwd=c(1,1),col=c("black","red"))

# create a plot where y-axis is percent increased coverage