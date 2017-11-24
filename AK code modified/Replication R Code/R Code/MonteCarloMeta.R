# Monte Carlo Meta Studies
rm(list=ls())
pathname<-"/Users/virasemenora/Dropbox (MIT)/Isaiah Andrews Replication"
setwd(pathname)
setwd("R Code")
library(xtable)
source("VariationVarianceLogLikelihood.R")
source("AuxiliaryFunctions.R")
# Number of MC simulations
nsims<-500



# True mean
lambdabar<-0;
# True variance
tau<-0.299;
# Number of draws in 1 replications
nbase<-10^4;
# Max data size
n<-100;
ndraws<-10^4;
critval<-1.96;
cutoffs=c(1.96);
symmetric<-1;
beta=c(0.15,1)
#Psi.true<-c(lambdabar,tau,beta[1])
alpha.true<-0.05
#pathname<-"/Users/virasemenora/Dropbox (MIT)/Isaiah Andrews Replication"
filepath<-"/Applications/EconExperiments/cleaned_econ_data.csv";
my.data<-read.csv(paste0(pathname,filepath),header=FALSE)
X<-my.data[,1]
sigma_econ<-my.data[,2]
sigma_scale<-mean(sigma_econ^2)

theta_hat_store<-matrix(0,nsims,2+(1-symmetric));
se_store<-matrix(0,nsims,2+(1-symmetric));
norm_theta<-matrix(0,nsims,2+(1-symmetric));
LM_store<-matrix(0,nsims,2+(1-symmetric));
if(FALSE) {
  theta_hat_store_add<-matrix(0,nsims,2+(1-symmetric));
  se_store_add<-matrix(0,nsims,2+(1-symmetric));
  norm_theta_add<-matrix(0,nsims,2+(1-symmetric));
  
}

Psihat0<-c(2,1)
Psi.true<-c(tau,beta[1])
for (ind in 1:nsims) {
  set.seed(ind)
  sigma<-rchisq(nbase, 2, ncp = 0)*sqrt(sigma_scale/8)
  
  # see Def. 3 Meta-study dgp, page 15
  Lambdastar<-lambdabar + tau*rnorm(nbase);
  X= Lambdastar+sigma*rnorm(nbase);
  TT=X/sigma;

  if (symmetric == 0) {
    Tpowers = cbind(TT<critval, TT>=critval);
  } else {
    Tpowers = cbind(abs(TT)<critval, abs(TT)>=critval);
  }
 
  D = (runif(nbase) < Tpowers %*%beta)
  
  # Trim data
  X = X[D];
  sigma = sigma[D];
  X= X[1:n]
  sigma = sigma[1:n]
  
 # thetahat0=c(lambdabar, tau, beta[1]);
  if (symmetric == 0) {
    LLH_only <-function (Psi) {
      A<-VariationVarianceLogLikelihood(Psi[1], Psi[2], c(Psi[-c(1,2)],  1),cutoffs,symmetric, X, sigma,0);
      return (A$LLH)
    }
    lower.b = c(-Inf,0,-Inf)
    upper.b=c(Inf,Inf,Inf) 
  } else {
    LLH_only <-function (Psi) {
      A<-VariationVarianceLogLikelihood(0, Psi[1], c(Psi[-c(1)],  1),cutoffs,symmetric, X, sigma,0);
      return (A$LLH)
    }
    lower.b = c(0,-Inf)
    upper.b=c(Inf,Inf)
  }
  
  findmin<-nlminb(objective=LLH_only, start=Psihat0,lower=lower.b,upper=upper.b);
  Psihat<-findmin$par
  LLH_max<-findmin$objective
  
  #findmin<-nlm(f=LLH_only,p=Psihat0,gradtol=0.00000001,
  #             iterlim=10000000)
  #Psihat<-findmin$estimate
  #LLH_max<-findmin$minimum
  stepsize = 0.001
  se = robust_se(Psihat,stepsize)
  theta_hat_store[ind,]<-Psihat
  se_store[ind,]<-se
  norm_theta[ind,] <- (Psihat-Psi.true)/se
  
  # remove restriction
}
maxval<-500
inds.na<-(1:maxval)[(is.na(se_store[,1]))]
theta_hat_store<-theta_hat_store[-inds.na,]
se_store<-se_store[-inds.na,]
norm_theta<-norm_theta[-inds.na,]

maxval<-dim(norm_theta)[1]
b.hat<-list()
#b.hat<-list(theta_hat_store=)
b.hat[["1"]]<-t(theta_hat_store[1:maxval,])
st.error.hat<-list()
st.error.hat[["1"]]<-t(se_store[1:maxval,])
norm.b.hat<-list()
norm.b.hat[["1"]]<-t(norm_theta[1:maxval,])





res<-list(b.hat = b.hat,
          st.error.hat = st.error.hat,
          norm.b.hat = norm.b.hat,
          which_coefs = 1:(2+1-symmetric))
filename<-"Results.tex"
directoryname<-paste0(pathname,"/Monte Carlo")
methods.my<-list("1")
#methods.my<-list("Meta")
names(methods.my)<-c("1")
add_info<-list(nsims,alpha.true,n)
names(add_info)<-c("Number of MC repetitions","Nominal test size", "Sample size")
output_result(res,
              true.value=Psi.true,filename=filename,directoryname,
              methods=methods.my,add_info=add_info,
              alpha = alpha.true,
              digs=3, N_rep=maxval)

