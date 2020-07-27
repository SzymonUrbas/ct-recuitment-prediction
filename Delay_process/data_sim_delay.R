## Simulate and add on delay
set.seed(1) # 1 or 100
C_true = 200
C = C_true
tauF = 900
tauH = tauF

if(F){
  tau_c = sort(c(tauF, sample(1:tauF,C-1, T)), T)
  
}
#tau_c = rep(tauF, C)
if(T){
  
  xx = c(1,rbeta(3*C/8,10,40),rbeta(C/8, 3,3),rbeta(C/4, 40,20), rbeta(C/8-1,45,5),1-rexp(C/8, 10))
  
  #tau_c = sort(ceiling(xx*tauF), T)
  
  #tau_c = sort(tauF-(rep(0:9/10*tauF, each = C/10)+rgeom(C, 0.1)),T)
  start_times_plan = sort(floor((1-xx)*600))
  tau_c = tauF-start_times_plan
  plot(tauF-tau_c,rep(0,C), pch = 3)
  
}
## quantiles of Weibull delay
p1 = 0.05; p2 = 0.95
q1 = 10; q2  =322
## solve for parameters
a_delay = log(log(1-p1)/log(1-p2))/log(q1/q2)
b_delay = q1/((-log(1-p1))^(1/a_delay))

## data:
#C = 60
delay_data = floor(rweibull(C,shape = a_delay, scale = b_delay))
start_times_true = start_times_plan+delay_data

idx_open = start_times_true<tauF
mean(idx_open)
C = sum(idx_open)
#start_times_plan = start_times_plan[idx]
start_times = start_times_true[idx_open]
tau_c = tauF - start_times
#tauF = floor(3*tauH/5)




start_times_x = start_times[start_times<tauF]
tau_c_x = tauF-start_times_x

tau_bar = mean(tau_c_x)

tauF = tauH
## true model - decaying
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
  
}



data = matrix(0, nrow = C, ncol = tauF)

N_cj = data

tilN_cj = data


for(c in 1:C){
  tilN_cj[c,1:(tau_c[c])] = rpois(tau_c[c],lambda0[c]*diff(c(0,G(1:(tau_c[c]),theta,tau_bar))))
  N_cj[c,1:(tau_c[c])] = tilN_cj[c,1:(tau_c[c])]
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
    tauF = 360#floor(3*tauH/5)
    #tauF = 200
    idx_x = as.logical(start_times<tauF)
    data_H = data[idx_x,]  # data available at the horizon
    data_F = data_H[,1:tauF]          # data available at the final observation day
    
    
    
    start_times_x = start_times[idx_x]
    tau_c_x = tauF-start_times_x
    
    idx_x = as.logical(c(idx_x,rep(0,C_true-C)))
    C_x = dim(data_F)[1]
    
    idx_censored = start_times_plan<tauF & (start_times_true)>=tauF
    
    N_cj_F = matrix(0, ncol = tauF, nrow = C_x)
    
    for(c in 1:C_x){
      
      N_cj_F[c,1:(tau_c_x[c])] = data_F[c,(start_times_x[c]+1):tauF]
    }
    
    tilN_cj = N_cj_F
    tilN_cj[,1] = N_cj_F[,1]-1 # remove first enrolment
    tilN_cj[(tilN_cj[,1]<0),1] = 0
  }
  delay_x = delay_data[idx_x]
  delay_censored = tauF - start_times_plan[idx_censored]
  
}
# par(mfrow = c(3,3))
# ind = sample(1:C, 9, replace = F)
# for(c in ind){
#  plot(cumsum(data[c,]), type = 'l', lwd = 2, ylim = c(0,20), xlab = 'Time', ylab = 'Accrual')
# }

par(mfrow = c(1,1))

