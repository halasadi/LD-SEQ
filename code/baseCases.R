#library(Hmisc) 
library(boot)

create_haps <- function(nsamp=100, nloci){
  # half are 1 0, and the other half are 0 1
  haps = matrix(nrow = nsamp, ncol = nloci, 0)
  for (i in 1:25){
    haps[i,1] = 1
    haps[i,2] = 0
  }
  for (i in 26:50){
    haps[i,1] = 0
    haps[i,2] = 0
  }
  for (i in 51:75){
    haps[i,1] = 1
    haps[i,2] = 1
  }
  for (i in 76:100){
    haps[i,1] = 0
    haps[i,2] = 1
  }
  return(haps)
}

pool <- function(lambda, haps, nloci, nsamp){
  ## simulate pool sequencing ##
  y_obs = rep(0, nloci)
  n = rep(0, nloci)
  n_1 = rep(0, nloci)
  for (i in 1:nloci){
    n[i] = rpois(1, lambda)
    n[i] = min(n[i], nsamp)
    sampled = sample(haps[,i], n)
    n_1[i] = sum(sampled)
    y_obs[i] = (n_1[i] + 0.5)/(length(sampled) + 1)
  }
  return(rbind(n, y_obs, n_1))
}

# try estimating sigma by ML

perform_estimate <- function(y_obs, nloci, nsamp, haps, pos){
  # calculate distance between SNPs
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
        dist = sum(d[min(i,j):(max(i,j)-1)])
        S[i,j] = exp(-rho*dist/nsamp)*cov_panel[i,j]
      }
    }
  }

  Sigma = ((1-theta)^2)*S + I((theta/2)*(1-(theta/2)))
  Sigma[abs(Sigma) < 1e-8] = 0


  d = diag(1/epsilon)
  Sigma_i = solve(Sigma) # there is also a problems, when nloci is big, kappa(Sigma) is very large (1e36)
  Sigma_bar = solve(Sigma_i + d)
  theta_bar = as.vector(Sigma_bar%*%(Sigma_i%*%mu + d%*%y_obs))
  return(as.vector(theta_bar))
}

sumsq <- function(x){
  return(sum(abs(x)))
}

nsamp = 100
nloci = 2
haps = create_haps(nsamp, nloci)
pos = c(50, 100)

nreps = 1000
lambdas = seq(5, 100, by = 5)
l = length(lambdas)

m_error_est = rep(0, l)
m_error_opt = rep(0, l)
m_error_obs  = rep(0, l)


for (ii in 1:l){
  temp_est = rep(0, nreps)
  temp_opt = rep(0, nreps)
  temp_obs = rep(0, nreps)
  for (jj in 1:nreps){
    lambda = lambdas[ii]
    pooled_info =  pool(lambda, haps, nloci, nsamp)
    n = pooled_info[1,]
    y_obs = pooled_info[2,]
    n_1 = pooled_info[3,]
    est = perform_estimate(y_obs, nloci, nsamp, haps, pos)
    true_freq = colMeans(haps)
    ## Matthew's suggestion is to do (n_1^1 + n_2^1)/(n_1+n_2), try it!
    opt_est = vector()
    opt_est[1] = (n_1[1] + (n[2]-n_1[2]))/ (n[1] + n[2])
    opt_est[2] = ((n[1]-n_1[1]) + n_1[2])/ (n[1] + n[2])
    temp_est[jj] = sumsq(est - true_freq)
    temp_opt[jj] = sumsq(opt_est - true_freq)
    temp_obs[jj] = sumsq(y_obs - true_freq)
  }
  m_error_est[ii] = mean(temp_est)
  m_error_opt[ii] = mean(temp_opt)
  m_error_obs[ii] = mean(temp_obs)

  # calculate variances here
}


plot(lambdas, m_error_est, xlab = "coverage", ylab = "sample mean of sum of LAE", col = "red", ylim = c(0,0.3))
points(lambdas, m_error_opt, col = "blue")
points(lambdas, m_error_obs)
legend("topright", c("read counts only at focal SNP","LDSP", "optimal bound"), lty=c(1,1,1), lwd=c(1,1,1),col=c("black","red", "blue"))
#legend("topright", c("read counts only at focal SNP","LDSP"), lty=c(1,1), lwd=c(1,1),col=c("black","red"))

## ratio of coverages
#plot(lambdas, m_error_est/m_error_obs, xlab = "true coverage", ylab = "effective coverage")

## next steps
## evolve and estimate sigma^2 by ML