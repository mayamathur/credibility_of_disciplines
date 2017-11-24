rm(list=ls())
pathname<-"/Users/virasemenora/Dropbox (MIT)/Isaiah Andrews Replication"
setwd(pathname)
setwd("R Code")
library(xtable)
source("ReplicationAnalyticLogLikelihood.R")
source("AuxiliaryFunctions.R")
# Number of MC simulations
nsims<-500


# True mean
lambdabar<-0;
# True variance
tau<-2.5;
# Number of draws in 1 replications
nbase<-10^4;
# Max data size
n<-50;
ndraws<-10^4;
critval<-1.96;
cutoffs=c(1.96);
symmetric=1;
betap=c(0.1,1);
alpha.true<-0.05



theta_hat_store<-matrix(0,nsims,2+(1-symmetric));
se_store<-matrix(0,nsims,2+(1-symmetric));
norm_theta<-matrix(0,nsims,2+(1-symmetric));
LM_store<-matrix(0,nsims,2+(1-symmetric));

for (ind in 1:nsims) {
  set.seed(ind+100)
Lambdastar<-lambdabar + tau*rnorm(nbase);
lambdasign=2*as.numeric(runif(nbase)<0.5)-1
#lambdasign<-1
Lambdastar=Lambdastar*lambdasign;
#normaldraws=matrix(rnorm(ndraws*2),ndraws,2)
#signdraws=as.numeric(runif(ndraws)<0.5)
Z = Lambdastar + rnorm(nbase);
# did not get this dividing by sigma 
#sigma = abs(Z)/3.242;
# page 9 says: use same density as for Z
sigma<-matrix(1,length(Z),1)
Zr = Lambdastar/sigma + rnorm(nbase);
if (symmetric == 1) {
  Tpowers = cbind(abs(Z)<critval, abs(Z)>=critval);
  
} else {
  Tpowers = cbind(Z<critval, Z>=critval);
  
}
D = (runif(nbase) < Tpowers %*%betap)

# Trim data
Z = cbind (Z[D], Zr[D]);
sigma = sigma[D];
Z = Z[1:n,];
sigma = sigma[1:n];
sigmaZ2 = sigma;


number=dim(Z)[1];
cluster_ID=1:number;
X=my.data[,1];
sigma<-my.data[,2];


#findmin<-nlm(f=LLH_only, p=Psihat0,steptol=0.00001,
#             gradtol = 0.00001,iterlim=100000000);

#Psihat<-findmin$estimate
#LLH_max<-findmin$minimum

if (symmetric==1) {
  Psihat0=c(2,1);
  Psi.true<-c(tau,betap[1])
  LLH_only <-function(Psi) {
   
    A<-ReplicationAnalyticLogLikelihood(0, Psi[1], c(Psi[-c(1)], 1),cutoffs,symmetric, Z, sigmaZ2);
    B<-A$LLH
    return(B)
  }
  lower.b = c(0,-Inf)
  upper.b=c(Inf,Inf)
  
} else {
  Psihat0=c(0,2,1);
  Psi.true<-c(lambdabar,tau,betap[1])
  LLH_only <-function(Psi) {
  #  Psi[3]<-betap[1]
 #   Psi[2]<-tau
    A<-ReplicationAnalyticLogLikelihood(Psi[1], Psi[2], c(Psi[-c(1,2)], 1),cutoffs,symmetric, Z, sigmaZ2);
    B<-A$LLH
    return(B)
  }
  lower.b = c(-Inf,0,-Inf)
  upper.b=c(Inf,Inf,Inf) 
}
findmin<-nlminb(objective=LLH_only, start=Psihat0,lower=lower.b,upper=upper.b);
Psihat<-findmin$par
LLH_max<-findmin$objective

theta_hat_store [ind,] = Psihat
stepsize = 0.001

se = robust_se(Psihat,stepsize);
se_store[ind,] = se
norm_theta[ind,] = (Psihat-Psi.true)/se


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

