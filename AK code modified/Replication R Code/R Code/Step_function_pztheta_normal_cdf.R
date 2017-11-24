#test<-FALSE
# tested both sym and asym
Step_function_pztheta_normal_cdf<-function (X,theta,sigma,alpha0, alpha1,cutoffs,symmetric) {
 # %arguments:
#    %X: point at which to evaluate cdf
 # %theta: parameter value under which to evaluate cdf
#  %sigma: standard deviation of (untruncated) normal variable
 # %sigma: stdev of distribution pi of mu
  #%alpha0: coefficient vector for intercept in logistic specification for p, in increasing order
  #%alpha1: coefficient vector for slope in logistic specification for p, in increasing order
  #%cutoffs:vector of thresholds for step function p, in increasing order.
  #%Cutoffs are given in terms of X, not z statistics
  #%symmetric: dummy indicating whether publication probability is symmetric
  #%around zero.  In symmetric case, cutoffs should include only positive
  #%values
  
  if (length(alpha0)!=(length(cutoffs)+1) || length(alpha1)!=(length(cutoffs)+1)) {
    stop('length of alpha0 and alpha1 must be one greater then length of cutoffs');
  }
  # %For symmetric case, create symmetrized version of cutoffs and coefficients
  
  if (symmetric ==1) {
    # cutoffs_u <-rep(0,2*length(cutoffs));
    # betap_u <-rep(0,2*length(cutoffs));
    
    cutoffs_u<-c(-rev(cutoffs),cutoffs);
    alpha0_u<-c(rev(alpha0)[1:length(cutoffs)],alpha0[1:length(cutoffs)]);
    alpha1_u<-c(rev(alpha1)[1:length(cutoffs)],alpha1[1:length(cutoffs)]);
    alpha0_u<-c(alpha0_u,1);
    alpha1_u<-c(alpha1_u,1);
  } else {
    cutoffs_u <-cutoffs;
    alpha1_u<-alpha1;
    alpha0_u<-alpha0;
    
  }
  
 # % %Calculate denominator in cdf
#  % denominator_calc=@(theta) sum(partition_probability_calculate(theta,1,cutoffs,symmetric)...
 #                                 %      .*(1+exp(-alpha0'*ones(1,length(theta'))-alpha1'*theta)).^-1,1);
                                  
  #                                %Calculate numerator in cdf
  
  cutoffs_u<-c(cutoffs_u,Inf);
  
  # %Calculate numerator in cdf
  if (X <= cutoffs_u[1]) {
    numerator=pnorm((X-theta)/sigma)/(1+exp(-alpha0_u[1]-theta*alpha1_u[1]));
  } else {
    numerator<-pnorm((cutoffs_u[1]-theta)/sigma)/(1+exp(-alpha0_u[1]-theta*alpha1_u[1]));
    m=1;
    while (X>cutoffs_u[m]) {
      Xcap=min(X,cutoffs_u[m+1]);
      
      numerator=numerator+(pnorm((Xcap-theta)/sigma)-pnorm((cutoffs_u[m]-theta)/sigma))/(1+exp(-alpha0_u[m]-theta*alpha1_u[m]));
      m<-m+1;
  }

  }
  
  
  Xtemp<-Inf;
 
  # Again?
cutoffs_u<-c(cutoffs_u,Inf);
  if (Xtemp<=cutoffs_u[1]) {
    denominator=pnorm((Xtemp-theta)/sigma)/(1+exp(-alpha0_u[1]-theta*alpha1_u[1]));
  
  } else {
    denominator=pnorm((cutoffs_u[1]-theta)/sigma)/(1+exp(-alpha0_u[1]-theta*alpha1_u[1]));
  }
   
  m<-1;
  while (Xtemp>cutoffs_u[m]) {
    Xcap=min(Xtemp,cutoffs_u[m+1]);
    denominator<-denominator+(pnorm((Xcap-theta)/sigma)-pnorm((cutoffs_u[m]-theta)/sigma))/(1+exp(-alpha0_u[m]-theta*alpha1_u[m]));
    m<-m+1;
  }
  
 

  cdf<-numerator/denominator;
  return(cdf)
}

#if (test) {
 # cutoffs=c(-1.96,1.96);
 
# alpha0 = c(-1.000, 0.0996,    1.0000);
 #alpha1 = c(-1.000,-0.5,    1.0000);
#  sigma=1;
#  theta=1;
#  X=1;
#  sigma<-1
#  symmetric<-0
#  cdf<-Step_function_pztheta_normal_cdf(1,theta,sigma,alpha0, alpha1,cutoffs,symmetric)
#}

#matcode= "/Users/virasemenora/Dropbox (MIT)/Isaiah Andrews Replication/Step_function_pztheta_normal_cdf.m"
