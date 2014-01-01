library(mvtnorm)

create_haps <- function(nsamp=100, nloci){
  haps = matrix(nrow = nsamp, ncol = nloci, 0)
  for (i in 1:25){
    haps[i,1] = 1
    haps[i,2] = 0
  }
  for (i in 26:50){
    haps[i,1] = 1
    haps[i,2] = 0
  }
  for (i in 51:60){
    haps[i,1] = 1
    haps[i,2] = 0
  }
  for (i in 61:100){
    haps[i,1] = 0
    haps[i,2] = 1
  }
  return(haps)
}

evolve_haps <- function(haps){
  # increase frequency of haplotype containing "1" allele at SNP 1
  # e.g. from 60 (1,0) haplotypes to 90 (1,0) haplotypes (a 50% increase in frequency)
  for (i in 61:90){
    haps[i,1] = 1
    haps[i,2] = 0
  }
  return(haps)
}


# to debug code
#global_obs = 0
#global_mu = 0
#global_Sigma = 0
#lobal_D = 0

findMLE <- function(obs, mu, Sigma, D){
  #global_obs <<- obs
  #global_mu <<- mu
  #global_Sigma <<- Sigma
  #global_D <<- D
  
  ll <- function(sigma2){
    cov_matrix = sigma2*Sigma + D
    sinv = solve(cov_matrix)
    p = length(mu)
    logl = -0.5*log(det(cov_matrix)) - 0.5*(t(obs-mu)%*%solve(cov_matrix)%*%(obs-mu)) -
      (p/2)*log(2*pi)
    # multiply by -1 since nlm minimizes
    return(as.numeric(-logl))
  }
  return(optim(ll, par = 1, method = "L-BFGS-B", lower = 0.1, upper = 10)$par)
}

pool <- function(lambda, haps, nloci, nsamp){
  ## simulate pool sequencing ##
  y_obs = rep(0, nloci)
  n = rep(0, nloci)
  n_1 = rep(0, nloci)
  for (i in 1:nloci){
    n[i] = rpois(1, lambda)
    # make sure it's below the number of sample in the pools
    n[i] = min(n[i], nsamp)
    # above 1
    n[i] = max(1, n[i])
    sampled = sample(haps[,i], n[i])
    n_1[i] = sum(sampled)
    y_obs[i] = (n_1[i] + 0.5)/(length(sampled) + 1)
  }
  return(rbind(n, y_obs, n_1))
}

# try estimating sigma by ML

perform_estimate <- function(y_obs, nloci, nsamp, haps, pos, n){
  # calculate distance between SNPs
  d = rep(0, nloci-1)
  for (i in 1:(nloci-1)){
    d[i] = as.numeric(pos[i+1])-as.numeric(pos[i])
  }

  # there is a problem here where n[i] = 0 for low coverage
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
        dist = sum(d[min(i,j):(max(i,j)-1)])
        S[i,j] = exp(-rho*dist/nsamp)*cov_panel[i,j]
      }
    }
  }

  Sigma = ((1-theta)^2)*S + I((theta/2)*(1-(theta/2)))
  Sigma[abs(Sigma) < 1e-8] = 0


  d = diag(1/epsilon)
  Sigma_i = solve(Sigma) # there is also a problems, when nloci is big, kappa(Sigma) is very large (1e36)
  
  sigma2 = findMLE(y_obs, mu, Sigma, diag(epsilon))
  
  # this is where the dispersion parameter comes in: eqn 13
  Sigma_bar = solve((1/sigma2)*Sigma_i + d)
  theta_bar = as.vector(Sigma_bar%*%((1/sigma2)*Sigma_i%*%mu + d%*%y_obs))
  return(as.vector(theta_bar))
}

mse <- function(y_hat, y){
  return(sum((y_hat-y)^2)/length(y))
}

nsamp = 100
nloci = 2
haps = create_haps(nsamp, nloci)
pos = c(50, 100)

nreps = 1000
lambdas = seq(5, 50, by = 5)
l = length(lambdas)

mse_est = rep(0, l)
mse_opt = rep(0, l)
mse_obs  = rep(0, l)

ev_haps = evolve_haps(haps)

for (ii in 1:l){
  store_ldsp_est = rep(0,nreps)
  store_opt_est = rep(0, nreps)
  store_obs_est = rep(0, nreps)
  store_true_freq = rep(0, nreps)
  for (jj in 1:nreps){
    lambda = lambdas[ii]
    pooled_info =  pool(lambda, ev_haps, nloci, nsamp)
    n = pooled_info[1,]  # total read counts
    n_1 = pooled_info[3,] # counts of "1" allele
    y_obs = pooled_info[2,] # estimate
    est = perform_estimate(y_obs, nloci, nsamp, haps, pos, n)
    true_freq = colMeans(ev_haps)
    store_true_freq[jj] = true_freq[1] # just the first SNP
    store_opt_est[jj] = (n_1[1] + (n[2]-n_1[2]))/ (n[1] + n[2])
    store_obs_est[jj] = y_obs[1]
    store_ldsp_est[jj] = est[1]
  }
  mse_est[ii] = mse(store_ldsp_est, store_true_freq)
  mse_opt[ii] = mse(store_opt_est, store_true_freq)
  mse_obs[ii] = mse(store_obs_est, store_true_freq)

}


plot(lambdas, mse_est, xlab = "coverage", ylab = "mean square error", col = "red", ylim = c(0,0.03), main = "frequency estimate for SNP 1", lwd=1.5)
points(lambdas, mse_opt, col = "blue", lwd=1.5)
points(lambdas, mse_obs, lwd=1.5)
legend("topright", c("read counts only at focal SNP","LDSP", "intuitive optimum"), lty=c(1,1,1), lwd=c(2,2,2),col=c("black","red", "blue"))
#legend("topright", c("read counts only at focal SNP","LDSP"), lty=c(1,1), lwd=c(1,1),col=c("black","red"))

# To-DO
# recheck code, make sure everything is solid!
# write evolve haps so you can specify percent increase
# make code more readable
# recheck code again
# investigate what the dispersion paramter should be optimally

## create a plot where y-axis is percent increased coverage