rm(list=ls())
set.seed(888)
source("ReplicationAnalyticLogLikelihood.R")
lambdabar<-0;
tau<-2.6297;
nsims<-500;
nbase<-10^3;
n<-18;
ndraws<-10^4;
critval<-1.96;
cutoffs=c(1.96);
symmetric=0;
beta=c(1,0.1);

theta_hat_store<-matrix(0,nsims,3);
se_store<-matrix(0,nsims,3);
se_robust_store<-matrix(0,nsims,3);
LM_store<-matrix(0,nsims,3);


for (k in 1:nsims) {
  set.seed(k)
 # %drawing latent parameters and simulation draws
  Lambdastar<-lambdabar + tau*rnorm(nbase);
  lambdasign=as.numeric(runif(nbase)<0.5)
  Lambdastar=Lambdastar*lambdasign;
  #normaldraws=matrix(rnorm(ndraws*2),ndraws,2)
  #signdraws=as.numeric(runif(ndraws)<0.5)
  Z = Lambdastar + rnorm(nbase);
  # did not get this dividing by sigma 
  #sigma = abs(Z)/3.242;
  # page 9 says: use same density as for Z
  sigma<-matrix(1,length(Z),1)
  Zr = Lambdastar/sigma + rnorm(nbase);
  Tpowers = cbind(abs(Z)>=critval, abs(Z)<critval);
  D = (rnorm(nbase) < Tpowers %*%beta)
  
  # Trim data
  Z = cbind (Z[D], Zr[D]);
  sigma = sigma[D];
  Z = Z[1:n,];
  sigma = sigma[1:n];
  sigmaZ2 = sigma;
  # Simulate Estimation
  # startvalue for estimation 
  # [lambdabar, tau, beta(2)]
  thetahat0 = c(0.02,2.5,0.3)
  
  #%create anonymous function that holds data and draws fixed
  LLH <-function(Psi) {
    A<-ReplicationAnalyticLogLikelihood(Psi[1], Psi[2], c(Psi[-c(1,2)], 1),cutoffs,symmetric, Z, sigmaZ2);
    A<-A$LLH
    return(A)
  }
  findmin<-nlm(f=LLH,p=thetahat0,gradtol = 0.0000001,steptol=0.000001,iterlim=1000000)
  thetahat<-findmin$par
  stepsize = 10^(-6)
  Info = matrix(0,length(thetahat),length(thetahat));
  for (n1 in 1:length(thetahat)) {
    for (n2 in 1:length(thetahat)) {
      thetaplusplus=thetahat;
      thetaplusminus=thetahat;
      thetaminusplus=thetahat;
      thetaminusminus=thetahat;
      
      thetaplusplus[n1]=thetaplusplus[n1]+stepsize;
      thetaplusplus[n2]=thetaplusplus[n2]+stepsize;
      LLH_plusplus=LLH(thetaplusplus);
      
      thetaplusminus[n1]=thetaplusminus[n1]+stepsize;
      thetaplusminus[n2]=thetaplusminus[n2]-stepsize;
      LLH_plusminus=LLH(thetaplusminus);
      
      thetaminusplus[n1]=thetaminusplus[n1]-stepsize;
      thetaminusplus[n2]=thetaminusplus[n2]+stepsize;
      LLH_minusplus=LLH(thetaminusplus);
      
      thetaminusminus[n1]=thetaminusminus[n1]-stepsize;
      thetaminusminus[n2]=thetaminusminus[n2]-stepsize;
      LLH_minusminus=LLH(thetaminusminus);
      
      
      Info[n1,n2]=((LLH_plusplus-LLH_plusminus)/(2*stepsize)-(LLH_minusplus-LLH_minusminus)/(2*stepsize))/(2*stepsize);
      
    }
  }
  Var = Info^(-1)
  se = sqrt(diag(Var))
  theta_hat_store [k,] = thetahat;
  se_store[k,] = se;
  #LLHmax<-findmin$
    
}
