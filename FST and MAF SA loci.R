N.e = 10^6
theta = 0.01 #4Neu
s.max = 0.01
u = theta/(4*N.e)

loci = 10^4


#simulations with drift
locus.counter = 0
F.ST = rep(0, loci)
MAF = rep(0, loci)
s.avg = rep(0, loci)

while(locus.counter < loci){
  s.f = runif(1)*s.max
  s.m = runif(1, s.f/(1 + s.f), s.f/(1 - s.f))
  
  p = rbeta(1, 4*N.e*u, 4*N.e*u)
  x = runif(1)
  
  if(min(p, 1-p) > 0.01){
    diff = s.f*p*(1 - p)/(2*(1 - (1 - p)*s.f)) + s.m*p*(1 - p)/(2*(1 - p*s.m))
    F.ST[locus.counter+1] = diff^2/(4*p*(1 - p))
    MAF[locus.counter+1] = min(p, 1 - p)
    s.avg[locus.counter+1] = ((s.f + s.m)/2)^2
    locus.counter = locus.counter + 1
  }else{locus.counter = locus.counter + 0}
}


par(mfrow = c(1, 3))
plot(MAF, F.ST, pch = 16, col = "grey", xlab = "MAF", ylab = "F.ST")
abline(lm(F.ST~MAF))

cov.drift = cov(MAF, F.ST)


#simulations with selection and drift
locus.counter = 0
F.ST = rep(0, loci)
MAF = rep(0, loci)
s.avg = rep(0, loci)

while(locus.counter < loci){
  s.f = runif(1)*s.max
  s.m = runif(1, s.f/(1 + s.f), s.f/(1 - s.f))
  p.eq = (s.f - s.m + s.f*s.m)/(2*s.f*s.m)

  p = rbeta(1, 4*N.e*u, 4*N.e*u)
  x = runif(1)
  
  if(x < exp(-N.e*s.m*s.f*(p.eq - p)^2) & min(p, 1-p) > 0.01){
    diff = s.f*p*(1 - p)/(2*(1 - (1 - p)*s.f)) + s.m*p*(1 - p)/(2*(1 - p*s.m))
    F.ST[locus.counter+1] = diff^2/(4*p*(1 - p))
    MAF[locus.counter+1] = min(p, 1 - p)
    s.avg[locus.counter+1] = ((s.f + s.m)/2)^2
    locus.counter = locus.counter + 1
  }else{locus.counter = locus.counter + 0}
}


plot(MAF, F.ST, pch = 16, col = "grey", xlab = "MAF", ylab = "F.ST")
abline(lm(F.ST~MAF))

cov.sel.drift = cov(MAF, F.ST)

#deterministic simulations
locus.counter = 0
F.ST = rep(0, loci)
MAF = rep(0, loci)
s.avg = rep(0, loci)

while(locus.counter < loci){
  s.f = runif(1)*s.max
  s.m = runif(1, s.f/(1 + s.f), s.f/(1 - s.f))
  p.eq = (s.f - s.m + s.f*s.m)/(2*s.f*s.m)
  
  p = p.eq
  
  if(min(p, 1-p) > 0.01){
    diff = s.f*p*(1 - p)/(2*(1 - (1 - p)*s.f)) + s.m*p*(1 - p)/(2*(1 - p*s.m))
    F.ST[locus.counter+1] = diff^2/(4*p*(1 - p))
    MAF[locus.counter+1] = min(p, 1 - p)
    s.avg[locus.counter+1] = ((s.f + s.m)/2)^2
    locus.counter = locus.counter + 1
  }else{locus.counter = locus.counter + 0}
}

plot(MAF, F.ST, pch = 16, col = "grey", xlab = "MAF", ylab = "F.ST")

cov.deterministic = cov(MAF, F.ST)
abline(lm(F.ST~MAF))

print(c(cov.drift,cov.sel.drift,cov.deterministic))
cov(s.avg, MAF)
