## Simulate
set.seed(100) # 1 or 100
C = 200
tauF = 600
tauH = tauF

if(F){
  tau_c = sort(c(tauF, sample(1:tauF,C-1, T)), T)
  
}

if(T){
  
  xx = c(1,rbeta(3*C/8,10,40),rbeta(C/8, 3,3),rbeta(C/4, 40,20), rbeta(C/8-1,45,5),1-rexp(C/8, 10))
  #xx = c(1,rbeta(3*C/8,40,10),rbeta(C/4, 30,20),rbeta(C/4, 20,10), rbeta(C/8-1,15,5))
  
  #hist(xx)
  tau_c = sort(ceiling(xx*tauF), T)
  
  #tau_c = sort(tauF-(rep(0:9/10*tauF, each = C/10)+rgeom(C, 0.1)),T)
  
  #plot(tauF-tau_c,rep(0,C), pch = 3)
  
}
start_times = tauF-tau_c
tauF = floor(3*tauH/5)


start_times_x = start_times[start_times<tauF]
tau_c_x = tauF-start_times_x

tau_bar = mean(tau_c_x)

tauF = tauH
## true model
if(T){
  
  theta = log(0.02)
  k = 2.7
  mod_sim = list(g = g_kappa_closure(k), G = G_kappa_closure(k),
               nllik = N_llik_closure(G_kappa_closure(k)),
               lprior = lprior_closure(k),nlpost = Nlpost_closure(k))
  G = mod_sim$G
  g = mod_sim$g
  tauF = 3/5*tauH
  ts = seq(0,tauH, length = 1e2)
  G_true = G(ts, theta, tau_bar)
  plot(ts, G_true, type = 'l' )
  
  g_true = g(ts, theta, tau_bar)
  plot(ts, g_true, type = 'l' )
  tauF=tauH
  if(T){
    alpha = 1.4
    phi = 0.01
    lambda0 = rgamma(C, alpha, rate = alpha/phi)
  }
  if(F){
    alpha1 = 1.4
    
    alpha2 = 1.4
    
    phi_bar = 0.01; r = 10
    phi1 = 2*phi_bar/(r+1)
    phi2 = r*phi1
    
    nu = 0.5
    lams = seq(0,0.2, length = 2e2)
    plot(lams,nu*dgamma(lams,alpha1,rate = alpha1/phi1)+(1-nu)*dgamma(lams,alpha2,rate = alpha2/phi2),
         type = 'l')
    cc = rbinom(1,C,nu)
    lambda0 = sample(c(rgamma(cc,alpha1,rate = alpha1/phi1),
                       rgamma(C-cc, alpha2, rate = alpha2/phi2)), C, F)
    
  }
}

if(F){
  alpha = 1.4
  phi = 0.01
  a = 1.5
  b = 30
  
  b*((a-1)/a)^(1/a)
  
  theta = c(a,b)
  tauF = 3/5*tauH
  G = function(t,theta,tauF){
    a = theta[1]; b = theta[2]
    val = tauF*pweibull(t,a,b)/pweibull(tauF,a,b)
    return(val)
  }
  
  g = function(t,theta,tauF){
    a = theta[1]; b = theta[2]
    val = tauF*dweibull(t,a,b)/pweibull(tauF,a,b)
    return(val)
  }
  
  G = Vectorize(G, vectorize.args = 't')
  
  g = Vectorize(g, vectorize.args = 't')
  
  ts = seq(0,tauH, length = 1e3)
  G_true = G(ts, theta, tau_bar)
  #plot(ts, G_true, type = 'l' )
  
  g_true = g(ts, theta, tau_bar)
  #plot(ts, g_true, type = 'l' )
  
  lambda0 = rgamma(C, alpha, rate = alpha/phi)
  tauF=tauH
  
}

if(F){
  alpha = exp(0.35)
  phi = exp(-4)
  a = 2.5
  b = 100
  tauF = 3/5*tauH
  L=32/15#/tauF
  G = function(t,theta,tauF){
    if(t<=tauF/8){
      return(t*L*3/4)
    }
    else if(t<=tauF*2/8){
      return(L*(t-tauF/8)+tauF/8*3/4*L)
    }
    else if(t<=tauF*4/8){
      return(L/2*(t-2*tauF/8) + tauF/8*(3/4*L+L))
    }
    else {
      return(L/4*(t-4*tauF/8) + tauF/8*(3/4*L+L+L))
    }
  }
  
  g = function(t,theta,tauF){
    
    if(t<=tauF/8){
      return(L*3/4)
    }
    else if(t<=tauF*2/8){
      return(L)
    }
    else if(t<=tauF*4/8){
      return(L/2)
    }
    else {
      return(L/4)
    }
  }
  
  G = Vectorize(G, vectorize.args = 't')
  
  g = Vectorize(g, vectorize.args = 't')
  
  ts = seq(0,tauH, length = 1e2)
  G_true = G(ts, theta, tau_bar)
  plot(ts, G_true, type = 'l' )
  ts = seq(0,tauH, length = 1e2)
  g_true = g(ts, theta, tau_bar)
  plot(ts, g_true, type = 'l', ylim = c(0,L) )
  
  lambda0 = rgamma(C, alpha, rate = alpha/phi)
  tauF=tauH
}

if(F){
  alpha = exp(0.35)
  phi = exp(-4)
  a = 2.5
  b = 100
  theta = c(a,b)
  tauF = 3/5*tauH
  G = function(t,theta,tauF){
    
    val = t
    return(val)
  }
  
  g = function(t,theta,tauF){
   
    val = 1
    return(val)
  }
  
  G = Vectorize(G, vectorize.args = 't')
  
  g = Vectorize(g, vectorize.args = 't')
  
  ts = seq(0,tauH, length = 1e2)
  G_true = G(ts, theta, tauF)
  plot(ts, G_true, type = 'l' )
  
  g_true = g(ts, theta, tauF)
  plot(ts, g_true, type = 'l' )
  
  lambda0 = rgamma(C, alpha, rate = alpha/phi)
  tauF=tauH
}


data = matrix(0, nrow = C, ncol = tauF)

N_cj = data




for(c in 1:C){
  N_cj[c,1:(tau_c[c])] = rpois(tau_c[c],lambda0[c]*diff(c(0,G(1:(tau_c[c]),theta,tau_bar))))
  N_cj[c,1] = N_cj[c,1]+1
  data[c,(tauF-tau_c[c]+1):tauF] = N_cj[c,1:(tau_c[c])]
}
accrual = cumsum(colSums(data))
plot(0:tauF, c(0,accrual), type = 'l', main = '', xlab = 'Time', ylab = 'Accrual', lwd = 2)
points(start_times,rep(0,C), pch = 3, col = 'red')
# plot(0,0, type = 'l', ylim = c(-1,max(N_cj)+1), xlim = c(0,tauF))
# for(c in 1:C){
#  points(1:(tau_c[c]),N_cj[c,1:(tau_c[c])])
# }
# 

if(T){
  tilN_cj = N_cj
  tauH = tauF
  data_H = data
  data_F = data
  tau_c_x = tau_c
  C_x = C
  if(T){
    #tauF = floor(1*tauH/4)
    tauF = floor(3*tauH/5)
    #tauF = 200
    
    data_H = data[start_times<tauF,]  # data available at the horizon
    data_F = data_H[,1:tauF]          # data available at the final observation day
    
    
    
    start_times_x = start_times[start_times<tauF]
    tau_c_x = tauF-start_times_x
    
    
    C_x = dim(data_F)[1]
    
    N_cj_F = matrix(0, ncol = tauF, nrow = C_x)
    
    for(c in 1:C_x){
      
      N_cj_F[c,1:(tau_c_x[c])] = data_F[c,(start_times_x[c]+1):tauF]
    }
    
    tilN_cj = N_cj_F
    tilN_cj[,1] = N_cj_F[,1]-1 # remove first enrolment
    tilN_cj[(tilN_cj[,1]<0),1] = 0
  }
  
}
# par(mfrow = c(3,3))
# ind = sample(1:C, 9, replace = F)
# for(c in ind){
#  plot(cumsum(data[c,]), type = 'l', lwd = 2, ylim = c(0,20), xlab = 'Time', ylab = 'Accrual')
# }

par(mfrow = c(1,1))

