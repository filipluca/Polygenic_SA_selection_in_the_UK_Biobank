N.e = 10^6
theta = 0.01 #4Neu
u = theta/(4*N.e)
k = 0.25 #gamma shape parameter
loci = 5*10^3

#simulations with strong selection
s.avg = 0.01
theta = s.avg/k
locus.counter = 0
F.ST = rep(0, loci)
MAF = rep(0, loci)

while(locus.counter < loci){
  t = rgamma(1, shape = k, scale = theta)
  
  p = rbeta(1, 4*N.e*u, 4*N.e*u)
  x = runif(1)
  
  if(x < exp(-N.e*t*p) & min(p, 1-p) > 0.01){
    diff = t*p*(1 - p)/(2*(1 - p*t))
    F.ST[locus.counter+1] = diff^2/(4*p*(1 - p))
    MAF[locus.counter+1] = min(p, 1 - p)
    locus.counter = locus.counter + 1
  }else{locus.counter = locus.counter + 0}
}

par(mfrow = c(1, 3))
plot(MAF, F.ST, pch = 16, col = "grey", xlab = "MAF", ylab = "F.ST", log = "x")
abline(lm(F.ST~MAF))
MAF.1 = MAF
F.ST.1 = F.ST
cov.1 = cov(MAF, F.ST)


#simulations with intermediate selection
s.avg = 0.001
theta = s.avg/k
locus.counter = 0
F.ST = rep(0, loci)
MAF = rep(0, loci)

while(locus.counter < loci){
  t = rgamma(1, shape = k, scale = theta)
  
  p = rbeta(1, 4*N.e*u, 4*N.e*u)
  x = runif(1)
  
  if(x < exp(-N.e*t*p) & min(p, 1-p) > 0.01){
    diff = t*p*(1 - p)/(2*(1 - p*t))
    F.ST[locus.counter+1] = diff^2/(4*p*(1 - p))
    MAF[locus.counter+1] = min(p, 1 - p)
    locus.counter = locus.counter + 1
  }else{locus.counter = locus.counter + 0}
}

plot(MAF, F.ST, pch = 16, col = "grey", xlab = "MAF", ylab = "F.ST", log = "x")
abline(lm(F.ST~MAF))
MAF.2 = MAF
F.ST.2 = F.ST
cov.2 = cov(MAF, F.ST)


#simulations with weak selection
s.avg = 0.0001
theta = s.avg/k
locus.counter = 0
F.ST = rep(0, loci)
MAF = rep(0, loci)

while(locus.counter < loci){
  t = rgamma(1, shape = k, scale = theta)
  
  p = rbeta(1, 4*N.e*u, 4*N.e*u)
  x = runif(1)
  
  if(x < exp(-N.e*t*p) & min(p, 1-p) > 0.01){
    diff = t*p*(1 - p)/(2*(1 - p*t))
    F.ST[locus.counter+1] = diff^2/(4*p*(1 - p))
    MAF[locus.counter+1] = min(p, 1 - p)
    locus.counter = locus.counter + 1
  }else{locus.counter = locus.counter + 0}
}

plot(MAF, F.ST, pch = 16, col = "grey", xlab = "MAF", ylab = "F.ST", log = "x")
abline(lm(F.ST~MAF))
MAF.3 = MAF
F.ST.3 = F.ST
cov.3 = cov(MAF, F.ST)

print(c(cov.1, cov.2, cov.3))

par(mfrow=c(1,3))
plot(MAF.1, F.ST.1, pch = 16, col = "grey", xlab = "MAF", ylab = "F.ST", ylim = c(0, 5*10^-11))
#abline(lm((F.ST.1)~MAF.1))
plot(MAF.2, F.ST.2, pch = 16, col = "grey", xlab = "MAF", ylab = "F.ST", ylim = c(0, 5*10^-11))
#abline(lm((F.ST.2)~MAF.2))
plot(MAF.3, F.ST.3, pch = 16, col = "grey", xlab = "MAF", ylab = "F.ST", ylim = c(0, 5*10^-11))
#abline(lm((F.ST.3)~MAF.3))

