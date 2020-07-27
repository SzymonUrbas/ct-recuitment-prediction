## Set up functions for model construction (TAU) - changing to tau_bar
## NOTE: unless stated otherwise, take "par" variable to be the logs of parameters
##
## Null homogeneous model (Anisimov)
if(T){
  # Note: The Anisimov model is parametrised
  #       to result in a nested model interpretation
  N_llik0 = function(par, N_cj, tau_c){
    tauF = mean(tau_c)
    C = length(tau_c)
    alpha = exp(par[1]); phi = exp(par[2])
    
    N_sj = colSums(N_cj)
    N_cs = rowSums(N_cj)
    
    val_1 = C*(alpha*log(alpha/phi)-lgamma(alpha))
    val_2 = sum(lgamma(alpha+N_cs) - (alpha+N_cs)*log(tau_c+alpha/phi) )
    
    ## Curve-only term (Note: \theta = 0)
    H = rep(1,max(tau_c))
    val_3 = sum(N_sj*log(H)) 
    
    
    return(-(val_1+val_2+val_3))
    
  }
  
  # Shape and integrated shape for completeness, \theta parameter doesn't do anything
  g0 = function(t,theta,tauF){
    if(t<0){return(0)}
    return(1)
  }
  g0 = Vectorize(g0,vectorize.args = 't')
  G0 = function(t,theta,tauF){
    if(t<0){return(0)}
    return(t)
  }
  G0 = Vectorize(G0,vectorize.args = 't')
  lprior0 = function(par){
    # ALPHA
    pi_a = dnorm(par[1], mean = 0.2, sd = 2, log = T)
    
    # PHI
    if(abs( par[2] ) < 8){
      pi_ph = -log(16)
    }
    else{
      pi_ph = -1e5
    }
    
    return(pi_a+pi_ph)
  }
  nlpost0 = function(par,N_cj,tau_c){
    val = N_llik0(par,N_cj,tau_c)-lprior0(par)
    if(is.na(val)){return(1e5)}
    return(val)
  }
}


#Inhomogeneous models 
if(T){
  # Intensity shape function of two parameters
  g_k = function(t,theta_til,k,tauF){
    if(t<0){return(0)}
    theta = exp(theta_til)
    return(tauF*(1+theta*t/k)^(-k)*(theta/k)*(1-k)/((1+theta*tauF/k)^(1-k)-1))
  }
  g_k = Vectorize(g_k, vectorize.args = 't')
    
  G_k = function(tau,theta_til,k, tauF){
    if(tau<0){return(0)}
    theta = exp(theta_til)
    return(tauF*((1+theta*tau/k)^(1-k)-1)/((1+theta*tauF/k)^(1-k)-1))
  }
  G_k = Vectorize(G_k, vectorize.args = 'tau')
  
  # Fixing the \kappa paramter
  g_kappa_closure = function(kappa){
    #g_kappa = function(t,theta_til,tauF)
    if(kappa == 1){
      g1 = function(t, theta_til, tauF){
        if(t<0){return(0)}
        theta = exp(theta_til)
        g = tauF*1/(1+theta*t)*theta/log(theta*tauF+1)
        return(g)
      }
      g1 = Vectorize(g1, vectorize.args = 't')
      return(g1)
    }
    else if(kappa == Inf){
      gInf = function(t, theta_til, tauF){
        if(t<0){return(0)}
        theta=exp(theta_til)
        g = tauF*exp(-theta*t)*theta/(1-exp(-theta*tauF))
        return(g)
      }
      gInf = Vectorize(gInf, vectorize.args = 't') 
      return(gInf)
    }
    else{
      g = function(t,theta_til,tauF){
        return(g_k(t,theta_til,kappa,tauF))
      }
      return(g)
    }
  }
  
  G_kappa_closure = function(kappa){
    if(kappa == 1){
      G1 = function(tau, theta_til, tauF){
        if(tau<0){return(0)}
        theta = exp(theta_til)
        G = tauF*log(theta*tau+1)/log(theta*tauF+1)
        return(G)
      }
      G1 = Vectorize(G1, vectorize.args = 'tau')
      return(G1)
    }
    else if(kappa == Inf){
      GInf = function(tau, theta_til, tauF){
        if(tau<0){return(0)}
        theta=exp(theta_til)
        G = tauF*(1-exp(-theta*tau))/(1-exp(-theta*tauF))
        return(G)
      }
      GInf = Vectorize(GInf, vectorize.args = 'tau')
      return(GInf)
    }
    else{
      G = function(t,theta_til,tauF){
        return(G_k(t,theta_til,kappa,tauF))
      }
      return(G)
    }
  }
  
  # Template for the NEGATIVE log-likelihood for a given G (up to a constant)
  N_llik_closure = function(G){
    N_llik = function(par, N_cj, tau_c){
      
      tauF = mean(tau_c) # largest centre open time
      C = length(tau_c) # number of centres
      alpha = exp(par[1]); phi = exp(par[2])
      THETA = par[3]
      
      N_sj = colSums(N_cj)
      N_cs = rowSums(N_cj)
      
      # alpha, phi term
      val_1 = C*(alpha*log(alpha/phi)-lgamma(alpha))
      
      # interaction term
      val_2 = sum(lgamma(alpha+N_cs) - (alpha+N_cs)*log(G(tau_c, THETA, tauF)+alpha/phi) )
      ## Curve-only term
      tau1 = 1:(max(tau_c)); tau0 = tau1-1
      H = G(tau1,THETA, tauF)-G(tau0, THETA, tauF)
      H[H==0] = 1e-100 # handles possible numerical instabilities
      
      val_3 = sum(N_sj*log(H)) 
      
      # putting it all together
      return(-(val_1+val_2+val_3))
    }
    return(N_llik)
  }
  # Template for the log-prior
  lprior_closure = function(kappa){
    lprior = function(par){
      # alpha (normal prior based on meta-analysis)
      pi_a = dnorm(par[1], mean = 0.2, sd = 2, log = T)
      
      # phi (diffuse prior)
      if(abs( par[2] ) < 8){
        pi_ph = -log(16)
      } else{
        pi_ph = -1e5
      }
      # theta (uniform prior put on the intensity ratio at time t0)
      
      if(T){
        t0 = 120 # larger t0 => vaguer prior
        if(kappa == Inf){
          pi_th = log(t0)+par[3]-t0*exp(par[3])
        }
        else{
          k = kappa
          pi_th = par[3]+log(t0)-(k+1)*log(1+t0/k*exp(par[3]))
        }
      }
      
      if(F){
        A=1.1;B=1.1
        t0 = 120 # larger t0 => vaguer prior
        if(kappa == Inf){
          pi_th = log(t0)+par[3]-t0*exp(par[3])+dbeta(exp(-t0*exp(par[3])),A,B, log = T)
        }
        else{
          k = kappa
          pi_th = par[3]+log(t0)-(k+1)*log(1+t0/k*exp(par[3]))+
            dbeta((1+t0/k*exp(par[3]))^(-k),A,B, log = T)
        }
      }
      
      return(pi_a+pi_ph+pi_th)
    }
  }
  # Template for NEGATIVE log-posterior
  Nlpost_closure = function(kappa){
    nlpost = function(par,N_cj,tau_c){
      val = (N_llik_closure(G_kappa_closure(kappa)))(par,N_cj,tau_c)-(lprior_closure(kappa))(par)
      if(is.na(val)){return(1e5)} # this handles numerical instability for large/small alpha/phi
      return(val)
    }
  }
}



# Build models as lists g_model
# g       := shape
# G       := integrated shape
# nllik   := negative log-likelihood
# lprior  := log-prior
# nlpost  := negative log-posterior
kappa = c(0.5,1,2,Inf)

mod_1 = list(g = g_kappa_closure(kappa[1]), G = G_kappa_closure(kappa[1]),
             nllik = N_llik_closure(G_kappa_closure(kappa[1])),
             lprior = lprior_closure(kappa[1]), nlpost = Nlpost_closure(kappa[1]))

mod_2 = list(g = g_kappa_closure(kappa[2]), G = G_kappa_closure(kappa[2]),
             nllik = N_llik_closure(G_kappa_closure(kappa[2])),
             lprior = lprior_closure(kappa[2]),nlpost = Nlpost_closure(kappa[2]))

mod_3 = list(g = g_kappa_closure(kappa[3]), G = G_kappa_closure(kappa[3]),
             nllik = N_llik_closure(G_kappa_closure(kappa[3])),
             lprior = lprior_closure(kappa[3]), nlpost = Nlpost_closure(kappa[3]))

mod_4 = list(g = g_kappa_closure(kappa[4]), G = G_kappa_closure(kappa[4]),
             nllik = N_llik_closure(G_kappa_closure(kappa[4])),
             lprior = lprior_closure(kappa[4]),nlpost = Nlpost_closure(kappa[4]))

mod_0 = list(g = g0, G = G0, nllik = N_llik0, lprior = lprior0,nlpost = nlpost0)

mod_list = list(mod_1,mod_2,mod_3,mod_4,mod_0)
# the null model is last by convention

