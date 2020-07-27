### MLE and model selection


# run mle for each model

Model_Fits = list()
NH_lliks = rep(0,length(mod_list)-1) # max log-liks of Non-Homogeneous models

## Maximum Likelihood Estimation (MLE)
for(iter in 1:(length(mod_list)-1)){
  par_init = rep(-1,3) ## NH models
  MLE = optim(par_init, mod_list[[iter]]$nllik, N_cj = tilN_cj, tau_c = tau_c_x, hessian = T)
  
  ## Checks in case of lack of identifiability (may occur for \kappa < 1)
  if (det(MLE$hessian)<=0.001|det(MLE$hessian)>=1e6) {
    NH_lliks[iter] = -1e6
  } else if (sqrt(diag(solve(MLE$hessian)))[3]>12) {
    NH_lliks[iter] = -1e6
  } else {
    NH_lliks[iter] = -MLE$value
  }
  
  
  Model_Fits[[iter]] = MLE
}
## NULL model (last in the list by convention)
if(T){
  par_init = rep(-1,2)
  MLE_0 = optim(par_init, mod_list[[iter+1]]$nllik, N_cj = tilN_cj, tau_c = tau_c_x, hessian = T)
  H_llik = -MLE_0$value
  Model_Fits[[iter+1]] = MLE_0
  
} 

## Likelihood ratio test and model selection
if(T){
  X_test = 2*(NH_lliks-H_llik) # X_test ~ \chi^2_1
 
  if(sum(X_test<=qchisq(0.95,1))==length(X_test)){
    best_model = length(mod_list) # If True, non of the models are significantly better than the null
  }
  else{
    best_model = max(which(X_test == max(X_test))) ## ties can occur for large \kappa (models give the same predictions then)
  }
  #best_model =  length(mod_list)
  mod_ML = Model_Fits[[best_model]]
  par_est = mod_ML$par
  if(length(par_est) == 2){
    par_est = c(par_est,0) # this handles the case where the best model is the NULL
  }
  Sig = solve(mod_ML$hessian)
  alpha = exp(par_est[1]); phi = exp(par_est[2]); theta = par_est[3]
  beta = alpha/phi
  
  G = mod_list[[best_model]]$G # fix best shape
  
}

#for(m in length(mod))

## Model validation
if(T){
  par(mfrow = c(1,2))
  ## Random Effect distribution (~ gamma)
  # quantities needed for random effect posteriors
  N_c = rowSums(tilN_cj)
  G_c = G(tau_c_x,theta,tauF)
  
  lambda0 = (alpha+N_c)/(beta+G_c)
  hist(lambda0[1:C_x],probability =  T, main = 'lambda^o distribution')
  lines(dgamma(seq(0,max(lambda0)),shape = alpha, rate = beta))
  
  
  ## Centre accrual distribution in first tauD days (~ neg-binomial)
  tauD = ceiling(tauF/2)
  inds = (tau_c_x >= tauD)
  N_true = rep(0,sum(inds))
  i = 1
  for(c in 1:C_x){
    if(inds[c]){
      N_true[i] = sum(tilN_cj[c,1:tauD])
      i = i+1
    }
  }
  N_theo = rnbinom(1e3, size = alpha, prob = alpha/(alpha+phi*G(tauD, par_est[3],tauF)))
  
  qqplot(N_theo, N_true, xlab = 'Theoretical', ylab = 'Observed',
         main = 'Recruitments in first half (QQ-plot)')
  lines(range(N_true), range(N_true))
}

## Distribution of accrual under the model (can be used as diagnostic)
if(T){
  GG = matrix(0,ncol = tauH, nrow = C) # matrix of integrated intensities of each centre at each day (entry is 0 is centre is not open)
  for(c in 1:C){
    GG[c,(start_times[c]+1):tauH] = diff(c(0,G(1:(tau_c)[c], theta, tauF)))
  }
  # Random effect posteriors
  N_c = c(rowSums(tilN_cj),rep(0, C-C_x))
  G_c = c(G(tau_c_x,theta,tauF),rep(0, C-C_x))
  
  alpha_post = alpha+N_c; beta_post = beta+G_c
  
  # Sampling from predictive
  if(T){
    M=1e3 # number of samples from the predictive
    
    pred_dist = matrix(0, nrow = M, ncol = tauH)
    # Samples from predictive for each day (bit of an overkill, also a bit slow)
    for(m in 1:M){
      lams = rgamma(C, shape = alpha_post, rate = beta_post)
      for(t in 1:tauH){
        
        pred_dist[m,t] = sum(rpois(C,lams*GG[,t]))
      }
    }
    
  }
  # Mean and Prediction Intervals of the accrual
  if(T){
    first_enrol = rep(0,tauH)
    for(ind in start_times){
      first_enrol[ind+1] = first_enrol[ind+1]+1 # there are probably neater ways of getting this
    }
    accrual_pred = t(apply(pred_dist, 1, cumsum))
    E_pred = apply(accrual_pred, 2, mean) + cumsum(first_enrol)
    
    pred_CI = apply(accrual_pred,2, quantile, probs = c(0.025,0.975))
    pred_CI[1,]=pred_CI[1,] + cumsum(first_enrol)
    pred_CI[2,]=pred_CI[2,] + cumsum(first_enrol)
  
   
    
    
  }
  ## Plotting
  if(T){
    par(mfrow = c(1,1))
    plot(0:tauH, c(0,cumsum(colSums(data))),lwd = 2, type = 'l', xlab = 'Time', ylab = 'Accrual')
    lines(0:tauH, c(0, E_pred), col = 'red', lwd = 2)
    for(i in 1:2){
      lines(0:tauH,c(0,pred_CI[i,]), col = 'red', lty = 2)
    }
    abline(v=tauF, lty = 2)
    points(start_times,rep(0,C), pch = 3, col = 'red')
  }
  
  
}



## Forecast from tauF to tauH
if(tauF<tauH){
  # Mean and Prediction Intervals of the accrual
  if(T){
    first_enrol = rep(0,tauH)
    for(ind in start_times){
      first_enrol[ind+1] = first_enrol[ind+1]+1 # there are probably neater ways of getting this
    }
    accrual_fcast = cbind(rep(0,M),t(apply(pred_dist[,(tauF+1):tauH], 1, cumsum)))
    E_fcast = apply(accrual_fcast, 2, mean) + cumsum(first_enrol)[(tauF):tauH] + sum(N_c)
    
    fcast_CI = apply(accrual_fcast,2, quantile, probs = c(0.025,0.975))
    fcast_CI[1,]=fcast_CI[1,] + cumsum(first_enrol)[(tauF):tauH]+ sum(N_c)
    fcast_CI[2,]=fcast_CI[2,] + cumsum(first_enrol)[(tauF):tauH]+ sum(N_c)
    
    MLE_E_fcast  = E_fcast
    MLE_fcast_CI = fcast_CI
    
    
  }
  ## Plotting
  if(T){
    par(mfrow = c(1,1))
    plot(0:tauH, c(0,cumsum(colSums(data))),lwd = 2, type = 'l', xlab = 'Time', ylab = 'Accrual')
    lines(tauF:tauH, E_fcast, col = 'red', lwd = 2)
    for(i in 1:2){
      lines(tauF:tauH,fcast_CI[i,], col = 'red', lty = 2)
    }
    points(start_times,rep(0,C), pch = 3, col = 'red')
    points(tauF, E_fcast[1], col = 'red', lwd = 4, pch  =1)
  }
}
