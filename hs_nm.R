### methods: 
#### 1. slice: stepping-out slice sampling dealing with tau and lambda(others: gibbs sampling)
#### 2. MwG: Metropolis within Gibbs sampling (MH:lambda & tau)
#### 3. gibbs: gibbs sampling using Wand mixture
#### 4. MSwG: Metropolis slice within Gibbs sampling (MH:lambda; MH: tau)
#### 5. slice2: Damien et al.(1999) slice sampling dealing with tau and lambda (others: gibbs sampling)
hs_nm <- function(Y,iter=5000,burn_in=1000,method=c("slice","MwG","gibbs","MSwG","slice2")){
  ## initial
  T <- length(Y)
  sigma2 <- 0.95*var(Y)
  tau <- 1
  lambda <- rep(1,T)
  beta <- Y 
  eta <- 1
  v <- rep(1,T)

  
  ## functions
  ### sigma
  sigma_gibbs <- function(T,beta,Y,lambda,tau){
    sigma2 <- 1/rgamma(1,T,rate=0.5*sum((Y-beta)^2)+0.5*sum(beta^2/lambda^2/tau^2))
    return(sigma2)
  }
  
  ### beta: blocked gibbs
  beta_gibbs <- function(Y,sigma2,tau,lambda){
    tl <- tau^2*lambda^2
    Sig <- sqrt(sigma2*tl/(1+tl))
    mu <- tl*Y/(1+tl)
    beta <- rnorm(T,mu,Sig)
    return(beta)
  }
  
  ### tau
  tau_slice <- function(T,beta,lambda,tau,sigma2){
    ## log transformation
    ff <- function(x){
      #f <- exp(-0.5*sum(beta^2/lambda^2/sigma2)*x)*x^(0.5*(T-1))/(1+x)
      f <- exp(-0.5*sum(beta^2/lambda^2/sigma2)/exp(2*x))*exp(x)^(-(T-1)*x)/(1+exp(2*x))
      return(f)
    }
    
    #eta <- 1/tau^2
    eta <- log(tau)
    eta_last <- eta
    f <- ff(eta)
    u1 <- runif(1,0,f)
    ## stepping out
    w <- 0.2
    m <- 5
    ### initial interval
    L <- eta - w*runif(1,0,1)
    #L <- max(L,0)
    R <- L + w
    ### limit of extension
    J <- floor(m*runif(1,0,1))
    K <- (m-1) - J
    while (J>0 & u1 < ff(L)){
      L <- L - w
      #L <- max(L,0)
      J <- J - 1
    }
    while (K>0 & u1 < ff(R)){
      R <- R + w
      K <- K - 1
    }
    ## shrinkage
    Lbar <- L
    Rbar <- R
    Accept <- 0
    while (!Accept){
      u2 <- runif(1,0,1)
      eta <- Lbar + u2*(Rbar-Lbar)
      if (u1 < f){
        Accept <- 1
      }else{
        if (eta < eta_last){
          Lbar <- eta
        }else{
          Rbar <- eta
        }
      }
    }
    #if (is.nan(sqrt(eta))){
    if (is.nan(exp(eta))){
      tau <- tau
      #}else{tau <- 1/sqrt(eta)}
    }else{tau <- exp(eta)}
    return(tau)
  }
  tau_slice2 <- function(T,beta,lambda,tau,sigma2){
    eta <-  1/tau^2
    u1 <- runif(1,0,1/(1+eta))
    u2 <- runif(1,0,pgamma(1/u1-1,(T+1)/2, 0.5*sum(beta^2/lambda^2/sigma2)))
    eta <- qgamma(u2,(T+1)/2, 0.5*sum(beta^2/lambda^2/sigma2))
    if (is.nan(sqrt(1/eta))){
      tau <- tau
    }else{
      tau <- 1/sqrt(eta)
    }
    return(tau)
  }
  tau_mh <- function(T,beta,lambda,tau,sigma2){
    ## inverse square transformation
    eta <- 1/tau^2
    y <- eta + rnorm(1,0,0.2)
    y <- abs(y)
    deta <- exp(-0.5*sum(beta^2/lambda^2/sigma2)*eta)*eta^(0.5*(T-1))/(1+eta)
    dy <- exp(-0.5*sum(beta^2/lambda^2/sigma2)*y)*y^(0.5*(T-1))/(1+y)
    alpha <- min(1,dy/deta)
    u <- runif(1,0,1)
    if (is.nan(sqrt(y))){
      tau <- tau
    }else{
      if(u<alpha){
        tau <- 1/sqrt(y)
      }else{
        tau <- tau
      }
    }
    return(tau)
  }
  tau_wm <- function(T,beta,lambda,sigma2,eta){
    tau2 <- 1/rgamma(1,0.5*(T+1),rate=1/eta+0.5*sum(beta^2/lambda^2)/sigma2)
    tau <- sqrt(tau2)
    return(tau)
  }

  ### lambda
  lambda_slice <- function(T,beta,tau,lambda,sigma2){
    ## inverse square transformation
    xi <- 1/lambda^2
    for (j in 1:T){
      f <- 1/(1+xi)
      mu2 <- beta[j]^2/tau^2/sigma2/2
      u1 <- runif(1,0,f)
      u2 <- runif(1,0,1-exp(-mu2*(1/u1-1)))
      xi[j] <- -log(1-u2)/mu2
      if (is.nan(sqrt(xi[j]))){
        lambda[j] <- lambda[j]
      }else{lambda[j] <- 1/sqrt(xi[j])}
    }
    return(lambda)
  }
  lambda_mh <- function(T,beta,tau,lambda,sigma2){
    dxi <- rep(0,T)
    dy1 <- rep(0,T)
    ## inverse square transformation
    xi <- 1/lambda^2
    alpha <- rep(0,T)
    y1 <- rep(0,T)
    for (j in 1:T){
      y1[j] <- xi[j] + rnorm(1,0,0.2)
      y1[j] <- abs(y1[j])
      dxi[j] <- exp(-0.5*beta[j]^2/tau^2/sigma2*xi[j])/(1+xi[j])
      dy1[j] <- exp(-0.5*beta[j]^2/tau^2/sigma2* y1[j])/(1+ y1[j])
      alpha[j] <- min(1,dy1[j]/dxi[j])
      u <- runif(1,0,1)
      if (is.nan(sqrt(y1[j]))){
        lambda[j] <- lambda[j]
      }else{
        if(u<alpha[j]){
          lambda[j] <- 1/sqrt(y1[j])
        }else{
          lambda[j] <- lambda[j]
        }
      }
    }
    return(lambda)
  }
  lambda_wm <- function(T,beta,tau,sigma2,v){
    lambda2 <- 1/rgamma(T,1,rate=1/v+0.5*(beta^2/tau^2)/sigma2)
    lambda <- sqrt(lambda2)
    return(lambda)
  }

  
  ## main
  sigma_sample <- rep(0,iter)
  beta_sample <- matrix(0,nrow=iter,ncol=T)
  tau_sample <- rep(0,iter)
  lambda_sample <- matrix(0,nrow=iter,ncol=T)
  for (i in 1:iter){
    ### update beta
    beta <- beta_gibbs(Y,sigma2,tau,lambda)
    
    ### update tau and lambda
    if (method=="slice"){
      tau <- tau_slice(T,beta,lambda,tau,sigma2)
      lambda <- lambda_slice(T,beta,tau,lambda,sigma2)
    }else if (method=="MwG"){
      tau <- tau_mh(T,beta,lambda,tau,sigma2)
      lambda <- lambda_mh(T,beta,tau,lambda,sigma2)
    }else if (method=="gibbs"){
      tau <- tau_wm(T,beta,lambda,sigma2,eta)
      eta <- 1/rgamma(1,1,rate=1+1/tau^2)
      lambda <- lambda_wm(T,beta,tau,sigma2,v)
      v <- 1/rgamma(T,1,rate=1+1/lambda^2)
    }else if (method=="MSwG"){
      tau <- tau_slice(T,beta,lambda,tau,sigma2)
      lambda <- lambda_mh(T,beta,tau,lambda,sigma2)
    }else if (method=="slice2"){
      tau <- tau_slice2(T,beta,lambda,tau,sigma2)
      lambda <- lambda_slice(T,beta,tau,lambda,sigma2)
    }
    
    ### update sigma2
    sigma2 <- sigma_gibbs(T,beta,Y,lambda,tau)
    
    sigma_sample[i] <- sigma2
    beta_sample[i,] <- beta
    tau_sample[i] <- tau
    lambda_sample[i,] <- lambda
  }
  
  ## burn_in
  sigma_sample <- sigma_sample[-(1:burn_in)]
  beta_sample <- beta_sample[-(1:burn_in),]
  tau_sample <- tau_sample[-(1:burn_in)]
  lambda_sample <- lambda_sample[-(1:burn_in),]
  
  
  
  ## return
  return(list(sigma2=sigma_sample,beta=beta_sample,tau=tau_sample,lambda=lambda_sample))
}



