## call ms ##
mutation_rate = 1e-6
popsize = 1e6
theta = 4*mutation_rate*popsize
nsamp = 100
r = 1e-8
nsites = 1e4
rho =  4 * popsize * r
mscall <- paste( "~/bin/msdir/ms", nsamp, 1, "-t", theta, "-r", rho*(nsites-1), nsites, "| sed '1,5d' > msout") 
system(mscall) # avoid this for relatively high recombination rate, will make R freeze

## read ms input ##
# read in positions
system("grep ^positions msout > pos.txt")
pos = read.table("pos.txt", header = FALSE, sep = "")
pos[1] = NULL
pos = floor(pos*nsites)
# read in haplotypes
raw_haps = system("cat msout | sed '1,1d'", intern=TRUE)

## format ms output ##
nloci = nchar(raw_haps[1])

haps = matrix(nrow = nsamp, ncol = nloci, 0)
for (i in 1:nsamp){
  str = raw_haps[i]
  for (j in 1:nloci){
    haps[i,j] = as.numeric(substring(str,j ,j))
  }
}


## simulate pool sequencing ##
lambda = 10  # coverage
y_obs = rep(0, nloci)
n = rep(0, nloci)

for (i in 1:nloci){
  n[i] = rpois(1, lambda)
  sampled = sample(haps[,i], n)
  y_obs[i] = mean(sampled)
}



# this part is tricky
# delete information about SNPs that have a pooled frequency of 0
y_temp = y_obs
pos_temp = pos
n_temp = n
haps_temp = haps

pos = n = y_obs = vector()
haps = NULL

cnt = 1
tol = 1e-6
for (i in 1:nloci){
  if (y_temp[i] > tol){
    y_obs[cnt] = y_temp[i]
    pos[cnt] = pos_temp[i]
    n[cnt] = n_temp[i]
    haps = cbind(haps, haps_temp[,i])
    cnt = cnt + 1
  }
}


# calculate distance between SNPs
nloci = length(pos)
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
Sigma[Sigma < 1e-8] = 0


d = diag(1/epsilon) # there is a problem here, often epsilons are 0
Sigma_i = solve(Sigma) # there is also a problems, when nloci is big, kappa(Sigma) is very large (1e36)

Sigma_bar = solve(Sigma_i + d)
theta_bar = Sigma_bar%*%(Sigma_i%*%mu + d%*%y_obs)
true_freq = colMeans(haps)
error_est = abs(true_freq-theta_bar)
error_obs = abs(true_freq-y_obs)
plot(error_est, error_obs, xlab = expression(paste("|", f_(true), "-", f_(estimate), "|")), ylab = expression(paste("|", f_(true), "-", y_obs, "|")))
lines(abs(true_freq-theta_bar), abs(true_freq-theta_bar))
plot(true_freq, error_est/error_obs, xlab = "true frequency", ylab ="error ratio (< 1 method does better)")
abline(h = mean(error_est/error_obs), col = "red")



