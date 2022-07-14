N.e = 10^6
mean.het = 0.01 #4Neu
u = mean.het/(4*N.e)
k = 1 #gamma shape parameter
r.fm = 0.9 #correlation coefficient
loci = 5*10^3

f.SL = 0.2 #fraction of loci that are sex-limited

#simulations with strong selection
s.avg = 0.01
theta = s.avg/k
locus.counter = 0
F.ST = rep(0, loci)
MAF = rep(0, loci)

while(locus.counter < loci){
  tf = rgamma(1, shape = k, scale = theta)
  y = rpois(1, tf*r.fm/(theta*(1 - r.fm)))
  tm = rgamma(1, shape = k + y, scale = theta*(1 - r.fm))
  
  tm = tm*rbinom(1, 1, (1 - f.SL))
  p = rbeta(1, 4*N.e*u, 4*N.e*u)
  x = runif(1)
  
  if(x < exp(-N.e*(tf + tm)*p) & min(p, 1-p) > 0.01){
    p.f = p + tf*p*(1 - p)/(2*(1 - p*tf))
    p.m = p + tm*p*(1 - p)/(2*(1 - p*tm))
    diff = p.f - p.m
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
  tf = rgamma(1, shape = k, scale = theta)
  y = rpois(1, tf*r.fm/(theta*(1 - r.fm)))
  tm = rgamma(1, shape = k + y, scale = theta*(1 - r.fm))
  
  tm = tm*rbinom(1, 1, (1 - f.SL))
  p = rbeta(1, 4*N.e*u, 4*N.e*u)
  x = runif(1)
  
  if(x < exp(-N.e*(tf + tm)*p) & min(p, 1-p) > 0.01){
    p.f = p + tf*p*(1 - p)/(2*(1 - p*tf))
    p.m = p + tm*p*(1 - p)/(2*(1 - p*tm))
    diff = p.f - p.m
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
  tf = rgamma(1, shape = k, scale = theta)
  y = rpois(1, tf*r.fm/(theta*(1 - r.fm)))
  tm = rgamma(1, shape = k + y, scale = theta*(1 - r.fm))
  
  tm = tm*rbinom(1, 1, (1 - f.SL))
  p = rbeta(1, 4*N.e*u, 4*N.e*u)
  x = runif(1)
  
  if(x < exp(-N.e*(tf + tm)*p) & min(p, 1-p) > 0.01){
    p.f = p + tf*p*(1 - p)/(2*(1 - p*tf))
    p.m = p + tm*p*(1 - p)/(2*(1 - p*tm))
    diff = p.f - p.m
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
