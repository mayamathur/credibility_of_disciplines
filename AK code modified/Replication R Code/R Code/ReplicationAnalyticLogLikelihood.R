# tested sym version
library(mvtnorm)
test<-FALSE
ReplicationAnalyticLogLikelihood <-function(nuhat, 
                                            tauhat, 
                                             betap, 
                                            cutoffs, 
                                            symmetric, 
                                            Z, 
                                            sigmaZ2) {
  
  ##show(nuhat)
  ##show(tauhat)
  ##show(betap)
  ##show(cutoffs)
  ##show(symmetric)
  ##show(Z)
  ##show(sigmaZ2)
  ##show("All!")

  
  browser()
  
  # MM note: I THINK NUHAT IS THE MEAN WE WANT TO CHANGE! 
  #arguments: mean and stdev of distribution pi of mu
  
  #coefficient vector for step function p, in increasing order
  #vector of thresholds for step function p, in increasing order
  #dummy for symmetric step function: if one, all cutoffs should be positive
  #n by 2 matrix of estimates
  #vector of stdevs of X2
  
  ###show(betap)
  n<-dim(Z)[1];
  k<-length(betap);
  betap<-as.matrix(betap,length(betap),1);
  #betap<-t(betap);
  Z1_dummies<-matrix(0,n,length(cutoffs)+1);
  
   #show(length(betap))
  if (all(sort(cutoffs)!=cutoffs)) {
    stop ("Unsorted cutoffs!")
  }
  
  if (symmetric ==1) {
    # check that cutoffs are positive
    if (!all(cutoffs>=0)) {
      stop("Needs positive cutoffs!")
    }
    Z1_dummies[,1]<-abs(Z[,1])<cutoffs[1];
    if (length(cutoffs)>1) {
      for (m in (2:length(cutoffs))) {
        Z1_dummies[,m]<-(abs(Z[,1])<cutoffs[m])*(abs(Z[,1])>=cutoffs[m-1]);
      }
    }
    Z1_dummies[,length(cutoffs)+1]<-abs(Z[,1])>=cutoffs[length(cutoffs)];
  } else {
    Z1_dummies[,1]=Z[,1]<cutoffs[1];
    if (length(cutoffs)>1) {
      for (m in 2:length(cutoffs)) {
        Z1_dummies[,m]=(Z[,1]<cutoffs[m])*(Z[,1]>=cutoffs[m-1]);
      }
    }
    Z1_dummies[,length(cutoffs)+1]<-Z[,1]>=cutoffs[length(cutoffs)];
    
  }

  
  #calculating components of the likelihood
 # #####show(dim(Z1_dummies))
  phat<-Z1_dummies%*%betap
  
  
  #% %likelihoods calculated by numerical integration
  #% piZ1Z2=zeros(n,1); %vector of un-truncated likelihoods
 # % for i=1:n
  #%     piZ1Z2(i)=2*normpdf(Z(i,1)-mudraws)'*1/sigmaZ2(i)*normpdf((Z(i,2)-mudraws)/sigmaZ2(i))/ndraws;
#% end
  # %likelihoods calculated analytically
  # %vector of un-truncated likelihoods
  piZ1Z2=rep(0,n); 
  for (i in 1:n) {
    Omega <- c( 1+tauhat^2, tauhat^2, tauhat^2, sigmaZ2[i]^2+tauhat^2 );
    Omega<-matrix(Omega,2,2);
    piZ1Z2[i]<-dmvnorm(c(Z[i,1], Z[i,2]),nuhat*c(1,1),Omega,log=FALSE)*0.5 +
      dmvnorm(c(Z[i,1], Z[i,2]),-nuhat*c(1,1),Omega,log=FALSE)*0.5
      
  }
 # %normalizing constant
  mu<-nuhat;
  sigma<-sqrt(tauhat^2+1);
 # %calculate probability that T statistic based on X_i has absolute value
#  %less than critval
  prob_vec<-rep(0,length(cutoffs)+2);
  
  if (symmetric==1) {
    for ( m in 1:length(cutoffs)) {
      prob_vec[m+1]<-0.5*(pnorm((cutoffs[m]-mu)/sigma)-pnorm((-cutoffs[m]-mu)/sigma))+
      0.5*(pnorm((cutoffs[m]+mu)/sigma)-pnorm((-cutoffs[m]+mu)/sigma));
    }
    prob_vec[length(cutoffs)+2]<-1;
    mean_Z1<-prob_vec[2:length(prob_vec)]-prob_vec[1:(length(prob_vec)-1)];
  } else {
    
    
    for ( m in 1:length(cutoffs)) {
    #prob_vec[m+1]<-0.5*(pnorm((cutoffs[m]-mu)/sigma)+pnorm((cutoffs[m]+mu)/sigma));
      prob_vec[m+1]<-0.5*(pnorm((cutoffs[m]-mu)/sigma)+pnorm((cutoffs[m]+mu)/sigma));
    }
    prob_vec[length(cutoffs)+2]<-1;
    mean_Z1<-prob_vec[2:length(prob_vec)]-prob_vec[1:(length(prob_vec)-1)];
  } 
  normalizingconst=mean_Z1%*%betap;

###show(betap)
  ###show(mean_Z1)
  ###show(piZ1Z2)
  ###show(phat)
  ###show(normalizingconst)
  ###show("All")
  #%vector of likelihoods
  # %cf equation 5 in the note
  L = as.vector(phat) * as.vector(piZ1Z2) / normalizingconst; 
  ####show(log(L))
  ####show(normalizingconst)
  if (normalizingconst<0) {
    ####show(betap)
  }
 
  logL=log(L);
  #  %objective function; note the sign flip, since we are doing minimization

  LLH = - sum(log(L));
  
  return(list(LLH=LLH,logL=logL))

}






### Testing 
if (test) {
  

nuhat<-5
tauhat<-1
betap<-c(0.1,0.2,0.3,0.4,0.5)
cutoffs<-c(-0.1503876, -0.2202119, 0.6473313, 0.9394543)
symmetric<-0
set.seed(1)
Z<-matrix(rnorm(2*100),100,2)
Z<-cbind(Test$V1,Test$V2)
sigmaZ2<-rep(1,100)


LLH<-ReplicationAnalyticLogLikelihood(nuhat, 
  tauhat, 
  betap, 
  cutoffs, 
  symmetric, 
  Z, 
  sigmaZ2)

}