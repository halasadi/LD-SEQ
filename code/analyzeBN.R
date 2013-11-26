nsamp = 100
theta = 5
mscall <- paste( "~/bin/msdir/ms", nsamp, 1, "-t", theta, "| sed '1,6d'") 
raw_haps = system( mscall, intern=TRUE)
nloci = nchar(raw_haps[1])
haps = matrix(nrow = nsamp, ncol = nloci, 0)

for (i in 1:nsamp){
  str = raw_haps[i]
  for (j in 1:nloci){
    haps[i,j] = as.numeric(substring(str,j ,j))
  }
}

