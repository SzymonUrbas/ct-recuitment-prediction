printf <- function(...){ invisible(print(sprintf(...)))}

## NH tests - individual

## simulate data:

## Likelihood-ratio =====================================================
# set s=c(1,1) when just 2 counts
pois.LRT =function(n.s, s){
  Fc = sum(s)
  n1s = n.s[1:s[1]]
  n2s = n.s[(s[1]+1):Fc]
  lam1 = sum(n1s)/s[1]
  lam2 = sum(n2s)/s[2]
  if(lam1 <= lam2){    return(0.5)  }
  if(lam2 == 0){return(0)}
  
  logLAM = sum(dpois(n.s, sum(n.s)/sum(s), log = T)) - sum(dpois(n1s,lam1,log = T))-sum(dpois(n2s,lam2,log = T))
  
  #logLAM = sum(n)*log(sum(n))-sum(n*log(n))+sum(n*log(s))
  
  p.val = 0.5*(1-pchisq(-2*logLAM, df = 1))
  return(as.numeric(p.val))
}  
LRT.power = function(lambda, s, alpha, iter = 1e3){
  p.val = rep(0,iter)
  for(i in 1:iter){
    # simulate data
    k = c(rpois(s[1], lambda[1]),rpois(s[2], lambda[2]))
    p.val[i] = pois.LRT(k,s)
  }
  return(list(power = sum(p.val<=alpha)/iter, dist = p.val))
}
LRT.power.table = function(lambdas, s, R, alpha, iter = 1e3){
  pows = matrix(0,length(lambdas), length(R))
  for(i in 1:length(lambdas)){
    for(j in 1:length(R)){
      pows[i,j] = LRT.power(lambda = lambdas[i]*c(1,R[j]), s = s, alpha = alpha, iter = iter)$power
      printf('Row: %2d out of %2d || Column: %2d out of %2d', i, length(lambdas), j, length(R))
    }
    
  }
  return(pows)
}


##
#==========================================================================
##



### ----------------------- non-par BS test ----------------------------- ###
### tests the difference


NPBS.test = function(n.s,s, BS = 1e4){
  
  obs.diff = sum(n.s[1:s[1]])/s[1]-sum(n.s[(s[1]+1):(sum(s))])/s[2]
  if(obs.diff<=0){return(1)}
  dist = replicate(BS,sum(sample(n.s,s[1],T))/s[1]-sum(sample(n.s,s[2],T))/s[2])
  return(sum(dist>=obs.diff)/BS)
}

NPBS.power = function(lambda, s, alpha, iter = 1e3, BS = 1e3){
  p.val = replicate(iter, NPBS.test(c(rpois(s[1],lambda[1]), rpois(s[2],lambda[2])),s, BS = BS))
  return( list( power = sum(p.val<=alpha)/iter, dist =  p.val))
}

NPBS.power.table  = function(lambdas, s, R, alpha, iter = 1e3){
  pows = matrix(0,length(lambdas), length(R))
  for(i in 1:length(lambdas)){
    for(j in 1:length(R)){
      pows[i,j] = NPBS.power(lambda = lambdas[i]*c(1,R[j]), s = s, alpha = alpha, iter = iter)$power
      printf('Row: %2d out of %2d || Column: %2d out of %2d', i, length(lambdas), j, length(R))
    }
  }
  return(pows)
}
