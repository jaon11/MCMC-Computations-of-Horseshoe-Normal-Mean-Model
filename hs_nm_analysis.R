#setwd()
library(ggplot2)
library(patchwork)
library(horseshoe)
source("hs_nm.R")
n <- 1
sigmean.out1 <- rep(0,n)
sigmean.out2 <- rep(0,n)
sigmean.out3 <- rep(0,n)
sigmean.out4 <- rep(0,n)
sigmean.out5 <- rep(0,n)
sigmean.out6 <- rep(0,n)
sigvar.out1 <- rep(0,n)
sigvar.out2 <- rep(0,n)
sigvar.out3 <- rep(0,n)
sigvar.out4 <- rep(0,n)
sigvar.out5 <- rep(0,n)
sigvar.out6 <- rep(0,n)

taumean.out1 <- rep(0,n)
taumean.out2 <- rep(0,n)
taumean.out3 <- rep(0,n)
taumean.out4 <- rep(0,n)
taumean.out5 <- rep(0,n)
taumean.out6 <- rep(0,n)
tauvar.out1 <- rep(0,n)
tauvar.out2 <- rep(0,n)
tauvar.out3 <- rep(0,n)
tauvar.out4 <- rep(0,n)
tauvar.out5 <- rep(0,n)
tauvar.out6 <- rep(0,n)

TP <- matrix(0,n,6)
FP <- matrix(0,n,6)
FN <- matrix(0,n,6)
TN <- matrix(0,n,6)


for (i in 1:n){
  set.seed(2003+i)
  signal <- 5
  mu <- c(rep(0,10),rep(signal,10),rep(0,10))
  #mu <- c(rep(0,10),rep(signal,10),rep(0,10),rep(signal/2,10),rep(0,10))
  eps <- rnorm(length(mu),0,1) # white noise
  # eps <- arima.sim(model = list(ar = 0.1), n = length(mu)) # AR(1) noise
  Y <- mu + eps
  T <- length(Y)
  out1 <- hs_nm(Y,method = "gibbs",iter = 5000,burn_in = 1000)
  out2 <- hs_nm(Y,method = "slice",iter = 5000,burn_in = 1000)
  out3 <- hs_nm(Y,method = "MwG",iter = 5000,burn_in = 1000)
  out5 <- hs_nm(Y,method = "MSwG",iter = 5000,burn_in = 1000)
  out6 <- hs_nm(Y,method = "slice2",iter = 5000,burn_in = 1000)
  
  ## horseshoe package
  out4 <- HS.normal.means(Y,method.tau = "halfCauchy", method.sigma = "Jeffreys")
 
  ## sigma2
  sigmean.out1[i] <- mean(out1$sigma2)
  sigmean.out2[i] <- mean(out2$sigma2)
  sigmean.out3[i] <- mean(out3$sigma2)
  sigmean.out4[i] <- out4$Sigma2Hat
  sigmean.out5[i] <- mean(out5$sigma2)
  sigmean.out6[i] <- mean(out6$sigma2)
  sigvar.out1[i] <- var(out1$sigma2)
  sigvar.out2[i] <- var(out2$sigma2)
  sigvar.out3[i] <- var(out3$sigma2)
  sigvar.out4[i] <- var(out4$Sigma2Samples)
  sigvar.out5[i] <- var(out5$sigma2)
  sigvar.out6[i] <- var(out6$sigma2)
  
  ## ACF
  acf.sig1 <- acf(out1$sigma2,plot = FALSE)
  acf.sig2 <- acf(out2$sigma2,plot = FALSE)
  acf.sig3 <- acf(out3$sigma2,plot = FALSE)
  acf.sig4 <- acf(out4$Sigma2Samples,plot = FALSE)
  acf.sig5 <- acf(out5$sigma2,plot = FALSE)
  acf.sig6 <- acf(out6$sigma2,plot = FALSE)
  
  sacfdf1 <- with(acf.sig1,data.frame(lag=lag,acf=acf))
  sacfdf2 <- with(acf.sig2,data.frame(lag=lag,acf=acf))
  sacfdf3 <- with(acf.sig3,data.frame(lag=lag,acf=acf))
  sacfdf4 <- with(acf.sig4,data.frame(lag=lag,acf=acf))
  sacfdf5 <- with(acf.sig5,data.frame(lag=lag,acf=acf))
  sacfdf6 <- with(acf.sig6,data.frame(lag=lag,acf=acf))
  
  s1 <- ggplot(sacfdf1,aes(x=lag,y=acf)) + 
    geom_segment(aes(xend=lag,yend=0)) +
    geom_hline(yintercept = 0.05,linetype = "dashed",color="blue") + 
    geom_hline(yintercept = -0.05,linetype = "dashed",color="blue") + 
    ggtitle("gibbs")
  s2 <- ggplot(sacfdf2,aes(x=lag,y=acf)) +
    geom_segment(aes(xend=lag,yend=0)) +
    geom_hline(yintercept = 0.05,linetype = "dashed",color="blue") + 
    geom_hline(yintercept = -0.05,linetype = "dashed",color="blue") + 
    ggtitle("slice")
  s3 <- ggplot(sacfdf3,aes(x=lag,y=acf)) +
    geom_segment(aes(xend=lag,yend=0)) +
    geom_hline(yintercept = 0.05,linetype = "dashed",color="blue") + 
    geom_hline(yintercept = -0.05,linetype = "dashed",color="blue") + 
    ggtitle("MwG")
  s4 <- ggplot(sacfdf4,aes(x=lag,y=acf)) +
    geom_segment(aes(xend=lag,yend=0)) +
    geom_hline(yintercept = 0.05,linetype = "dashed",color="blue") + 
    geom_hline(yintercept = -0.05,linetype = "dashed",color="blue") + 
    ggtitle("horseshoe")
  s5 <- ggplot(sacfdf5,aes(x=lag,y=acf)) +
    geom_segment(aes(xend=lag,yend=0)) +
    geom_hline(yintercept = 0.05,linetype = "dashed",color="blue") + 
    geom_hline(yintercept = -0.05,linetype = "dashed",color="blue") + 
    ggtitle("MSwG")
  s6 <- ggplot(sacfdf6,aes(x=lag,y=acf)) +
    geom_segment(aes(xend=lag,yend=0)) +
    geom_hline(yintercept = 0.05,linetype = "dashed",color="blue") + 
    geom_hline(yintercept = -0.05,linetype = "dashed",color="blue") + 
    ggtitle("slice2")
  s <- s1 + s2 + s6 + s3 + s4 + s5
  s + plot_annotation(
    title = "ACF of sigma2",
    caption = "blue dashed lines are 95% confidence interval"
  )
  ggsave(s,filename = paste0("ACF_sigma2_",i,".png"),width = 10,height = 10)
  
  
  
 
  ## tau
  taumean.out1[i] <- mean(out1$tau)
  taumean.out2[i] <- mean(out2$tau)
  taumean.out3[i] <- mean(out3$tau)
  taumean.out4[i] <- out4$TauHat
  taumean.out5[i] <- mean(out5$tau)
  taumean.out6[i] <- mean(out6$tau)
  tauvar.out1[i] <- var(out1$tau)
  tauvar.out2[i] <- var(out2$tau)
  tauvar.out3[i] <- var(out3$tau)
  tauvar.out4[i] <- var(out4$TauSamples)
  tauvar.out5[i] <- var(out5$tau)
  tauvar.out6[i] <- var(out6$tau)

  
  ## ACF
  acf.tau1 <- acf(out1$tau,plot = FALSE)
  acf.tau2 <- acf(out2$tau,plot = FALSE)
  acf.tau3 <- acf(out3$tau,plot = FALSE)
  acf.tau4 <- acf(out4$TauSamples,plot = FALSE)
  acf.tau5 <- acf(out5$tau,plot = FALSE)
  acf.tau6 <- acf(out6$tau,plot = FALSE)
  
  
  tacfdf1 <- with(acf.tau1, data.frame(lag, acf))
  tacfdf2 <- with(acf.tau2, data.frame(lag, acf))
  tacfdf3 <- with(acf.tau3, data.frame(lag, acf))
  tacfdf4 <- with(acf.tau4, data.frame(lag, acf))
  tacfdf5 <- with(acf.tau5, data.frame(lag, acf))
  tacfdf6 <- with(acf.tau6, data.frame(lag, acf))
  
  t1 <- ggplot(tacfdf1, aes(x=lag, y=acf)) + 
    geom_segment(aes(xend=lag, yend=0)) +
    geom_hline(yintercept=0.05, linetype="dashed",color="blue") + 
    geom_hline(yintercept=-0.05, linetype="dashed",color="blue") + 
    ggtitle("gibbs")
  t2 <- ggplot(tacfdf2, aes(x=lag, y=acf)) +
    geom_segment(aes(xend=lag, yend=0)) +
    geom_hline(yintercept=0.05, linetype="dashed",color="blue") + 
    geom_hline(yintercept=-0.05, linetype="dashed",color="blue") + 
    ggtitle("slice")
  t3 <- ggplot(tacfdf3, aes(x=lag, y=acf)) +
    geom_segment(aes(xend=lag, yend=0)) +
    geom_hline(yintercept=0.05, linetype="dashed",color="blue") + 
    geom_hline(yintercept=-0.05, linetype="dashed",color="blue") + 
    ggtitle("MwG")
  t4 <- ggplot(tacfdf4, aes(x=lag, y=acf)) +
    geom_segment(aes(xend=lag, yend=0)) +
    geom_hline(yintercept=0.05, linetype="dashed",color="blue") + 
    geom_hline(yintercept=-0.05, linetype="dashed",color="blue") + 
    ggtitle("horseshoe")
  t5 <- ggplot(tacfdf5, aes(x=lag, y=acf)) +
    geom_segment(aes(xend=lag, yend=0)) +
    geom_hline(yintercept=0.05, linetype="dashed",color="blue") + 
    geom_hline(yintercept=-0.05, linetype="dashed",color="blue") + 
    ggtitle("MSwG")
  t6 <- ggplot(tacfdf6, aes(x=lag, y=acf)) +
    geom_segment(aes(xend=lag, yend=0)) +
    geom_hline(yintercept=0.05, linetype="dashed",color="blue") + 
    geom_hline(yintercept=-0.05, linetype="dashed",color="blue") + 
    ggtitle("slice2")
  t <- t1 + t2 + t6 + t3 + t4 + t5 
  t + plot_annotation(
    title = "ACF of tau",
    caption = "The blue dashed lines represent the 95% confidence interval"
  )
  ggsave(t,filename=paste0("ACF_tau_",i,".png"),width=10,height=10)
  
  
  ## beta
  betamean.out1 <- colMeans(out1$beta)
  betamean.out2 <- colMeans(out2$beta)
  betamean.out3 <- colMeans(out3$beta)
  betamean.out4 <- out4$BetaHat
  betamean.out5 <- colMeans(out5$beta)
  betamean.out6 <- colMeans(out6$beta)
  betalb.out1 <- apply(out1$beta,2,quantile,0.025)
  betalb.out2 <- apply(out2$beta,2,quantile,0.025)
  betalb.out3 <- apply(out3$beta,2,quantile,0.025)
  betalb.out4 <- out4$LeftCI
  betalb.out5 <- apply(out5$beta,2,quantile,0.025)
  betalb.out6 <- apply(out6$beta,2,quantile,0.025)
  betaub.out1 <- apply(out1$beta,2,quantile,0.975)
  betaub.out2 <- apply(out2$beta,2,quantile,0.975)
  betaub.out3 <- apply(out3$beta,2,quantile,0.975)
  betaub.out4 <- out4$RightCI
  betaub.out5 <- apply(out5$beta,2,quantile,0.975)
  betaub.out6 <- apply(out6$beta,2,quantile,0.975)
  

  
  df.beta.out1 <- data.frame(mean=betamean.out1,lb=betalb.out1,ub=betaub.out1,mu=mu,Y=Y)
  df.beta.out2 <- data.frame(mean=betamean.out2,lb=betalb.out2,ub=betaub.out2,mu=mu,Y=Y)
  df.beta.out3 <- data.frame(mean=betamean.out3,lb=betalb.out3,ub=betaub.out3,mu=mu,Y=Y)
  df.beta.out4 <- data.frame(mean=betamean.out4,lb=betalb.out4,ub=betaub.out4,mu=mu,Y=Y)
  df.beta.out5 <- data.frame(mean=betamean.out5,lb=betalb.out5,ub=betaub.out5,mu=mu,Y=Y)
  df.beta.out6 <- data.frame(mean=betamean.out6,lb=betalb.out6,ub=betaub.out6,mu=mu,Y=Y)
  
  
  p1 <- ggplot(df.beta.out1,aes(x=1:T,y=mean)) + 
    geom_point(color="red") + 
    geom_point(aes(y=mu)) +
    geom_point(aes(y=Y),color="blue") +
    geom_ribbon(aes(ymin=lb,ymax=ub),alpha=0.3) + 
    ggtitle("gibbs")
  p2 <- ggplot(df.beta.out2,aes(x=1:T,y=mean)) +
    geom_point(color="red") + 
    geom_point(aes(y=mu)) +
    geom_point(aes(y=Y),color="blue") +
    geom_ribbon(aes(ymin=lb,ymax=ub),alpha=0.3) + 
    ggtitle("slice")
  p3 <- ggplot(df.beta.out3,aes(x=1:T,y=mean)) +
    geom_point(color="red") + 
    geom_point(aes(y=mu)) +
    geom_point(aes(y=Y),color="blue") +
    geom_ribbon(aes(ymin=lb,ymax=ub),alpha=0.3) + 
    ggtitle("MwG")
  p4 <- ggplot(df.beta.out4,aes(x=1:T,y=mean)) +
    geom_point(color="red") + 
    geom_point(aes(y=mu)) +
    geom_point(aes(y=Y),color="blue") +
    geom_ribbon(aes(ymin=lb,ymax=ub),alpha=0.3) + 
    ggtitle("horseshoe")
  p5 <- ggplot(df.beta.out5,aes(x=1:T,y=mean)) +
    geom_point(color="red") + 
    geom_point(aes(y=mu)) +
    geom_point(aes(y=Y),color="blue") +
    geom_ribbon(aes(ymin=lb,ymax=ub),alpha=0.3) + 
    ggtitle("MSwG")
  p6 <- ggplot(df.beta.out6,aes(x=1:T,y=mean)) +
    geom_point(color="red") + 
    geom_point(aes(y=mu)) +
    geom_point(aes(y=Y),color="blue") +
    geom_ribbon(aes(ymin=lb,ymax=ub),alpha=0.3) + 
    ggtitle("slice2")
  p <- p1 + p2 + p6 + p3 + p4 + p5
  p + plot_annotation(
    title = "Figure of Shrinkage Estimation",
    subtitle = "Black: true value, Red: estimated value, Blue: observed value"
  )
  ggsave(p,filename = paste0("hs_nm_",i,".png"),width=10,height=10)
  
 
  ## signal-noise detection
  ### out1
  TP[i,1] <- sum(mu!=0 & betalb.out1*betaub.out1>0)
  FP[i,1] <- sum(mu==0 & betalb.out1*betaub.out1>0)
  TN[i,1] <- sum(mu!=0 & betalb.out1*betaub.out1<=0)
  FN[i,1] <- sum(mu==0 & betalb.out1*betaub.out1<=0)
  ### out2
  TP[i,2] <- sum(mu!=0 & betalb.out2*betaub.out2>0)
  FP[i,2] <- sum(mu==0 & betalb.out2*betaub.out2>0)
  TN[i,2] <- sum(mu!=0 & betalb.out2*betaub.out2<=0)
  FN[i,2] <- sum(mu==0 & betalb.out2*betaub.out2<=0)
  ### out3
  TP[i,3] <- sum(mu!=0 & betalb.out3*betaub.out3>0)
  FP[i,3] <- sum(mu==0 & betalb.out3*betaub.out3>0)
  TN[i,3] <- sum(mu!=0 & betalb.out3*betaub.out3<=0)
  FN[i,3] <- sum(mu==0 & betalb.out3*betaub.out3<=0)
  ### out4
  TP[i,4] <- sum(mu!=0 & betalb.out4*betaub.out4>0)
  FP[i,4] <- sum(mu==0 & betalb.out4*betaub.out4>0)
  TN[i,4] <- sum(mu!=0 & betalb.out4*betaub.out4<=0)
  FN[i,4] <- sum(mu==0 & betalb.out4*betaub.out4<=0)
  ### out5
  TP[i,5] <- sum(mu!=0 & betalb.out5*betaub.out5>0)
  FP[i,5] <- sum(mu==0 & betalb.out5*betaub.out5>0)
  TN[i,5] <- sum(mu!=0 & betalb.out5*betaub.out5<=0)
  FN[i,5] <- sum(mu==0 & betalb.out5*betaub.out5<=0)
  ### out6
  TP[i,6] <- sum(mu!=0 & betalb.out6*betaub.out6>0)
  FP[i,6] <- sum(mu==0 & betalb.out6*betaub.out6>0)
  TN[i,6] <- sum(mu!=0 & betalb.out6*betaub.out6<=0)
  FN[i,6] <- sum(mu==0 & betalb.out6*betaub.out6<=0)
  
  
  
  
}
    
  

## computation effieciency
df.sig <- data.frame(mean=c(mean(sigmean.out1),mean(sigmean.out2),mean(sigmean.out6),mean(sigmean.out3),mean(sigmean.out4),mean(sigmean.out5)),
                        var=c(mean(sigvar.out1),mean(sigvar.out2),mean(sigvar.out6),mean(sigvar.out3),mean(sigvar.out4),mean(sigvar.out5)),
                        method=c("gibbs","slice","slice2","MwG","horseshoe","MSwG"))
print(df.sig)

## computation stability
df.tau <- data.frame(mean=c(mean(taumean.out1),mean(taumean.out2),mean(taumean.out6),mean(taumean.out3),mean(taumean.out4),mean(taumean.out5)),
                        var=c(mean(tauvar.out1),mean(tauvar.out2),mean(tauvar.out6),mean(tauvar.out3),mean(tauvar.out4),mean(tauvar.out5)),
                        method=c("gibbs","slice","slice2","MwG","horseshoe","MSwG"))
                               
print(df.tau)

## signal-noise detection
TP <- colSums(TP)
FP <- colSums(FP)
TN <- colSums(TN)
FN <- colSums(FN)

precision <- TP/(TP+FP)
recall <- TP/(TP+FN)
F1 <- 2*precision*recall/(precision+recall)
df <- data.frame(precision,recall,F1,method=c("gibbs","slice","MwG","horseshoe","MSwG","slice2"))
print(df)



