## Bayes model averaging (importance sampling) WITH DELAY
library(mvtnorm)

M_prop = 1e4 # number of proposal samples for each posterior
M      = 1e2 # number of samples from each posterior  (M <= M_prop)

Model_prior = c(1,1,1,1,1) # prior probabilities of each model


Post_samples = array(0,dim = c(M, 3, length(mod_list)))
# Dim = number of samples, number of parameters, number of models

tau_bar = mean(tau_c_x)

w_bar = rep(0,length(mod_list))
ESSs  = rep(0,length(mod_list))            

C_vec = rep(0,M)
# Get posterior modes and samples through importance sampling
# with resampling (makes averaging easier)

lw_max = rep(0,length(mod_list))


## ESTIMATE DELAY DISTRIBUTION
if(T){
  nllik_delay = function(par){
    a = exp(par[1]);b=exp(par[2])
    val = -sum(log(pweibull(delay_x+1,shape = a, scale = b)-
                     pweibull(delay_x,shape = a, scale = b)))
    val = val-sum(pweibull(delay_censored,shape = a, scale = b, lower.tail = FALSE, log.p = TRUE))
    return(val)
  }
  delay_par = exp(optim(c(0,5), nllik_delay, hessian = T)$par)
}


# Non-homogeneous models
for(iter in 1:(length(mod_list)-1)){
  # Posterior mode
  par_init = c(-1,-1,-2)
  print(mod_list[[iter]]$nllik(par_init,  N_cj = tilN_cj, tau_c = tau_c_x))
  print(mod_list[[iter]]$lprior(par_init))
  print(mod_list[[iter]]$nlpost(par_init,  N_cj = tilN_cj, tau_c = tau_c_x))
  p_mode = optim(par_init, mod_list[[iter]]$nlpost,  N_cj = tilN_cj, tau_c = tau_c_x, hessian = T)
  md = p_mode$par
  V = as.matrix(solve(p_mode$hessian))
  
  
  
  # t-proposal on 4 degrees of freedom
  
  prop_samp = rmvt(M_prop, sigma = V, delta = md, df = 4)
  
  # weighting ( log-weights )
  lw = rep(0, M_prop)
  for(i in 1:M_prop){
    lw[i] = -mod_list[[iter]]$nlpost(prop_samp[i,], tilN_cj, tau_c_x)-
      dmvt(prop_samp[i,], sigma = V, delta = md, df = 4, log = T) 
  }
  
  lw_max[iter] = max(lw)
  lw = lw-max(lw) # stabilises before exponentiating
  
  w = exp(lw)
  idx = sample(1:M_prop, M, prob = w, replace = T)
  Post_samples[,,iter] = prop_samp[idx,]
  w_bar[iter] = mean(w)
  ESSs[iter]  = 1 / sum( (w/sum(w))^2 )
  
}

# Null model (last in the list by convention)
if(T){
  # Posterior mode
  par_init = rep(0,2)
  print(mod_list[[iter+1]]$nlpost(par_init,  N_cj = tilN_cj, tau_c = tau_c_x))
  p_mode = optim(par_init, mod_list[[iter+1]]$nlpost,  N_cj = tilN_cj, tau_c = tau_c_x, hessian = T)
  md = p_mode$par
  V = as.matrix(solve(p_mode$hessian))
  
  
  # t-proposal on 4 degrees of freedom
  
  prop_samp = rmvt(M_prop, sigma = V, delta = md, df = 4)
  
  # weighting ( log-weights )
  lw = rep(0, M_prop)
  for(i in 1:M_prop){
    lw[i] = -mod_list[[iter+1]]$nlpost(prop_samp[i,], tilN_cj, tau_c_x)-
      dmvt(prop_samp[i,], sigma = V, delta = md, df = 4, log = T) 
  }
  
  lw_max[iter+1] = max(lw)
  lw = lw-max(lw) # stabilises before exponentiating
  
  w = exp(lw)
  idx = sample(1:M_prop, M, prob = w, replace = T)
  Post_samples[,1:2,iter+1] = prop_samp[idx,]
  w_bar[iter+1] = mean(w)
  ESSs[iter+1]  = 1 / sum( (w/sum(w))^2 )
}

# Calculating relative posterior model probabilities
log_EVIDs = log(w_bar)+lw_max
EVIDs = exp(log_EVIDs-max(log_EVIDs))
EVIDs = EVIDs/sum(EVIDs)
#EVIDs = c(0,0,0,0,1)

## model diagnostics are the same as in the MLE file


## Predictive distribution marginalised over all the models (including deviance)
if(T){
  mod_idx = sample(1:length(mod_list),M,prob = EVIDs*Model_prior, replace = T) # M "samples of models"
  
  N_c = c(rowSums(tilN_cj),rep(0, C-C_x))
  pred_dist = matrix(0, nrow = M, ncol = tauH)
  accrual_fcast = matrix(0, nrow = M, ncol = tauH-tauF+1)
  DEVs = matrix(0, nrow = M, ncol = 2)
  
  EN_t = pred_dist
  VN_t = pred_dist
  
  
  for(m in 1:M){
    if(m %% 10 == 0){
      print(m)
    }
    start_times_m = rep(0,C_true)
    start_times_m[idx_x] = start_times_x
    
    delay_m = floor(rweibull(C_true-C_x-sum(idx_censored), delay_par[1],  delay_par[2]))
    start_times_m[!(idx_x+idx_censored)] = start_times_plan[!(idx_x+idx_censored)]+delay_m 
    
    for( c in which(idx_censored)){
      delay = -1
      while(delay<tauF-start_times_plan[c]){
        delay = floor(rweibull(1,delay_par[1],delay_par[2]))
      }
      start_times_m[c] = start_times_plan[c]+delay
    }
    idx_pred = start_times_m<tauH
    C = sum(idx_pred)
    C_vec[m] = C
    start_times_m = start_times_m[idx_pred]
    tau_c = tauH-start_times_m
    
    preddata = matrix(0,nrow = C, ncol = tauH)
    predN_cj = matrix(0,nrow = C, ncol = tauH)
    

    first_enrol = rep(0,tauH)
    for(ind in start_times_m){
      first_enrol[ind+1] = first_enrol[ind+1]+1 # there are probably neater ways of getting this
    }
    # homogeneous model
    if(mod_idx[m]==length(mod_list)){
      mod_m = mod_list[[(mod_idx[m])]]
      ## insert stuff here
      alpha_m = exp(Post_samples[m,1,(mod_idx[m])])
      beta_m = alpha_m/exp(Post_samples[m,2,(mod_idx[m])])
      G_c = c(mod_m$G(tau_c_x,0,tau_bar),rep(0, C-C_x))
      
      GG = matrix(0, nrow = C, ncol = tauH)
      
      
      for(c in 1:C){
        GG[c,(start_times_m[c]+1):tauH] = diff(c(0,mod_m$G(1:(tau_c)[c], 0, tau_bar)))
      }
      
      
      beta_post = beta_m + G_c ; alpha_post = alpha_m + N_c
      lams = rgamma(C, shape = alpha_post, rate = beta_post)
      for(t in 1:tauH){
        p_t = beta_post/(beta_post+rowSums(as.matrix(GG[,1:t])))
        
        #EN_t[m,t] = sum( alpha_post*(1-p_t)/p_t )
        #VN_t[m,t] = sum( alpha_post*(1-p_t)^2/p_t )
        #sample rates lambda
        
        preddata[,t] = rpois(C,lambda = lams*GG[,t])
        #preddata[,t] = rnbinom(C,size = alpha_post,
        #                      prob = p_t)
        pred_dist[m,t] = sum(preddata[,t])
        
        
      }
      
      for(c in 1:C_x){
        predN_cj[c,1:(tau_c_x[c])] = preddata[c,(tauF-tau_c_x[c]+1):tauF]
      }
      
      DEVs[m,] = -c(mod_m$nllik(Post_samples[m,,(mod_idx[m])], predN_cj[1:C_x,1:tauF],tau_c_x),
                    mod_m$nllik(Post_samples[m,,(mod_idx[m])], tilN_cj,tau_c_x))
      
      
    }
    
    ## Inhomogeneous model
    else{
      mod_m = mod_list[[(mod_idx[m])]]
      
      alpha_m = exp(Post_samples[m,1,(mod_idx[m])])
      beta_m = alpha_m/exp(Post_samples[m,2,(mod_idx[m])])
      theta_m = Post_samples[m,3,(mod_idx[m])]
      
      G_c = c(mod_m$G(tau_c_x,theta_m,tau_bar),rep(0, C-C_x))
      
      GG = matrix(0, nrow = C, ncol = tauH)
      
      for(c in 1:C){
        GG[c,(start_times_m[c]+1):tauH] = diff(c(0,mod_m$G(1:(tau_c[c]), theta_m, tau_bar)))
      }
      
      
      beta_post = beta_m + G_c ; alpha_post = alpha_m + N_c
      lams = rgamma(C, shape = alpha_post, rate = beta_post)
      for(t in 1:tauF){
        
        # p_t = beta_post/(beta_post+rowSums(as.matrix(GG[,1:t])))
        # 
        # EN_t[m,t] = sum( alpha_post*(1-p_t)/p_t )
        # VN_t[m,t] = sum( alpha_post*(1-p_t)/(p_t^2) )
        preddata[,t] = rpois(C,lambda = lams*GG[,t])
        # preddata[,t] = rnbinom(C,size = alpha_post,
        # prob = p_t)
        pred_dist[m,t] = sum(preddata[,t])
        
      }
      for(t in (tauF+1):tauH){
        
        p_t = beta_post/(beta_post+rowSums(as.matrix(GG[,(tauF+1):t])))
        
        #EN_t[m,t] = sum( alpha_post*(1-p_t)/p_t )
        #VN_t[m,t] = sum( alpha_post*(1-p_t)/(p_t^2) )
        preddata[,t] = rpois(C,lambda = lams*GG[,t])
        # preddata[,t] = rnbinom(C,size = alpha_post,
        # prob = p_t)
        pred_dist[m,t] = sum(preddata[,t])
        
      }
      for(c in 1:C_x){
        predN_cj[c,1:(tau_c_x[c])] = preddata[c,(tauF-tau_c_x[c]+1):tauF]
      }
      DEVs[m,] = c(mod_m$nllik(Post_samples[m,,(mod_idx[m])], predN_cj[1:C_x,1:tauF],tau_c_x),
                   mod_m$nllik(Post_samples[m,,(mod_idx[m])], tilN_cj,tau_c_x))
      
      
      
      
    }
    accrual_fcast[m,] = cumsum(pred_dist[m,tauF:tauH])+cumsum(first_enrol)[tauF:tauH]+sum(tilN_cj)
    
    
  }
  
  
}


accrual_fcast_delay = accrual_fcast

## Forecasting only
if(T){
  if(tauF<tauH){
    # Mean and Prediction Intervals of the accrual
    if(T){
      
      #accrual_fcast = accrual_pred[,tauF:tauH]
      E_fcast = apply(accrual_fcast, 2, mean) 
      
      fcast_CI = apply(accrual_fcast,2, quantile, probs = c(0.025,0.975))
      
      BMA_E_fcast  = E_fcast
      BMA_fcast_CI = fcast_CI
    }
    ## Plotting
    if(T){
      par(mfrow = c(1,1))
      plot(0:tauH, c(0,cumsum(colSums(data))),lwd = 2, type = 'l', xlab = 'Time', ylab = 'Accrual',
           ylim = c(0,max(fcast_CI)), main = 'Fitted Weibull delay')
      lines(tauF:tauH, E_fcast, col = 'red', lwd = 2)
      for(i in 1:2){
        lines(tauF:tauH,fcast_CI[i,], col = 'red', lty = 2)
      }
      points(start_times_true,rep(0,C_true), pch = 3, col = 'red')
      points(tauF, E_fcast[1], col = 'red', lwd = 4, pch  =1)
    }
  }
}

## EXTRA

if(F){
  par(mfrow = c(2,1))
  plot(cumsum(first_enrol), type = 'l') # "day 1" enrolment
  plot(apply(pred_dist == 0,2,mean), type = 'l') # probability of making no recruitment on a given day
}

