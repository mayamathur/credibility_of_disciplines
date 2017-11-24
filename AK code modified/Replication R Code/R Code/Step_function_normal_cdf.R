# tested both sym and asym vers
#test<-FALSE

# MM: I added vectorization wrt theta because integrate function passes multiple thetas
# definitely changes answer given by integrate function
Step_function_normal_cdf = Vectorize( function(X,theta,sigma,betap,cutoffs,symmetric) {
  
  
 # %arguments:
#    %X: point at which to evaluate cdf (Z-SCORE)
#  %theta: parameter value under which to evaluate cdf
#  %sigma: standard deviation of (untruncated) normal variable
#  %sigma: stdev of distribution pi of mu
#  %betap: coefficient vector for step function p, in increasing order
#  %cutoffs:vector of thresholds for step function p, in increasing order.
 # %Cutoffs are given in terms of X, not z statistics  (MM: I THINK THIS IS WRONG GIVEN CODE IN HORIZONTALBARS.R)
 # %symmetric: dummy indicating whether publication probability is symmetric
 # %around zero.  In symmetric case, cutoffs should include only positive
  #%values
 # %NOTE: publication probability for largest category (i.e. for those points beyond largest cutoff) normalized to one.  
  if (length(betap) != length(cutoffs)+1) {
    stop('length of betap must be one greater than length of cutoffs');
  }
 # %For symmetric case, create symmetrized version of cutoffs and coefficients
  
  if (symmetric ==1) {
   # cutoffs_u <-rep(0,2*length(cutoffs));
   # betap_u <-rep(0,2*length(cutoffs));
    
    cutoffs_u<-c(-rev(cutoffs),cutoffs);
    betap_u<-c(rev(betap)[1:length(cutoffs)],betap[1:length(cutoffs)]);
    betap_u<-c(betap_u,1);
  } else {
    cutoffs_u <-cutoffs;
    betap_u<-betap;
    #betap_u<-c(betap_u,1);
  }
  # %Calculate denominator in cdf

  prob_vec = rep(0,length(cutoffs_u)+1);
  
  # MM: prob of being below each cutoff value based on theta and sigma
  # so we're assuming the Z-scores are N(theta, sigma)
 # ~~~~~~~ MAYBE WE HAVE TO CHANGE TAU IN MY INTEGRANDS BECAUSE THETA IS APPARENTLY THE TRUE Z-SCORES?
  for (m in 1:length(cutoffs_u)) {
    prob_vec[m+1] = pnorm( (cutoffs_u[m]-theta ) / sigma );
  }
 
  prob_vec<-c(prob_vec,1);

  # MM: prob of being between each pair of cutoffs, I think
  mean_Z1 <- prob_vec[2:length(prob_vec)] - prob_vec[1:(length(prob_vec)-1)];
 
  # MM: I think this is 1 / expectation
  # sum of terms for each cutoff value
  denominator = mean_Z1 %*% betap_u;
  
  cutoffs_u<-c(cutoffs_u,Inf);
 
 
  if (X <= cutoffs_u[1]) {
    # MM: numerator is the integrand in the CDF, I think
    # where betap_u is the p(x^tilde) term
   numerator <- pnorm((X-theta)/sigma)*betap_u[1];
   
  } else {
    numerator <- pnorm((cutoffs_u[1]-theta)/sigma)*betap_u[1];
    m=1;
    
    
    # MM: this is the integral in the CDF (pg. 21)
    # evaluated as a sum over all of the cutoffs
    # takes probability that X falls between each pair of cutoffs (diff in pnorms)
    #  and multiplies it by the relevant publication prob (betap_u) for that interval
    #  and adds these over all intervals
    while (X>cutoffs_u[m]) {
      Xcap<-min(X,cutoffs_u[m+1]);
      numerator <- numerator + ( pnorm((Xcap-theta)/sigma) - pnorm((cutoffs_u[m]-theta)/sigma) )*betap_u[m+1];
      m<-m+1;
    #  if (m > length(cut))
    }
   
   
  }
  
  cdf = numerator/denominator;
  return(cdf)
  
} , vectorize.args = "theta" )

### Testing 
# if (test) {
# cutoffs=c(-1.96,1.96);
# betap = c(0.0996,    1.0000);
# betap = c(-1,0.0996,    1.0000);
# sigma=1;
# theta=1;
# X=1;
# sigma<-1
# symmetric<-0
# cdf<-Step_function_normal_cdf(1,theta,sigma,betap,cutoffs,symmetric)
# }
# 
# 








