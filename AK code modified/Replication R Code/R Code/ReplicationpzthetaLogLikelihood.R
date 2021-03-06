# tested sym version that is the only one possible
# not tested
ReplicationpzthetaLogLikelihood<-function(nuhat, tauhat, betap, gammap, cutoffs,symmetric, 
                                          Z, sigmaZ2, vcutoff, Zdraws, Thetadraws,zdrawsd,thetadrawsd) {
  
 # %arguments: mean and stdev of distribution pi of mu
#  %coefficient vector beta for specification of p, in increasing order
 # %coefficient vector gamma for specification of p, in increasing order
#  %vector of thresholds for step function p, in increasing order
 # %dummy for symmetric step function: if one, all cutoffs should be positive
#  %n by 2 matrix of estimates
 # %vector of stdevs of X2
#  %cutoff for latent variable V in specification of p
 # %simulation draws for z
#  %simulation draws for theta
 # %standard deviation of z simulation draws
  #%standard deviation of theta simulation draws
  
  browser()
  
  # sanity check
 if (is.null(dim(Zdraws)) || is.null(dim(Thetadraws)) ) {
   Zdraws<-as.matrix(Zdraws)
   Thetadraws<-as.matrix(Thetadraws)
 }
  
  n<-length(Z[,1]);
  nz<-length(Zdraws[,1]);
  ntheta<-length(Thetadraws[,1]);
  k<-length(betap);
  betap<-as.matrix(betap,length(betap),1);
  gammap<-as.matrix(gammap,length(gammap),1);
  
  
  #%calculate p(z,theta) with elementwise operations
  
  
  if (symmetric == 1) {
    
    if (k==2 ) {
      pztheta<- function (z,theta) 
        {
       
        A<-(abs(z)<cutoffs[1])*(betap[1]-(pnorm(vcutoff-theta)-pnorm(-vcutoff-theta))*gammap[1])+(abs(z)>=cutoffs[1])*(betap[2]-(pnorm(vcutoff-theta)-pnorm(-vcutoff-theta))*gammap[2]);
        
        
      #  ##show(betap[1]-(pnorm(vcutoff-theta)-pnorm(-vcutoff-theta))*gammap[1])
      #  ##show((abs(z)<cutoffs[1])*(betap[1]-(pnorm(vcutoff-theta)-pnorm(-vcutoff-theta))*gammap[1]))
        return(A)
      }
    } else if (k==3) {
      
      pztheta=function (z,theta)  { 
        A<- (abs(z)<cutoffs[1])*(betap[1]-(pnorm(vcutoff-theta)-pnorm(-vcutoff-theta))*gammap[1])+(abs(z)>=cutoffs[1] & abs(z)<cutoffs[2])*(betap[2]-(pnorm(vcutoff-theta)-pnorm(-vcutoff-theta))*gammap[2])+(abs(z)>=cutoffs[2])*(betap[3]-(pnorm(vcutoff-theta)-pnorm(-vcutoff-theta))*gammap[3]);
        return(A)
      }
    } else {
      stop("code only specified for 1 or 2 cutoffs!")
    }
  } else {
    stop("code only specified for symmetric case")
  }
  
  Zdrawmat<-matrix(rep(Zdraws[,1],ntheta),ncol=ntheta);
  Thetadrawmat<-matrix(rep(Thetadraws[,1],nz),nrow=nz,byrow=TRUE);
  ##show(pztheta(Zdrawmat,Thetadrawmat))
  denominator<-mean(mean(pztheta(Zdrawmat,Thetadrawmat)
                         *exp(-0.5*(1-zdrawsd^-2)*Zdrawmat^2+Zdrawmat*Thetadrawmat-0.5*Thetadrawmat^2 
                              -0.5*(tauhat^-2-thetadrawsd^-2)*Thetadrawmat^2+Thetadrawmat*nuhat*tauhat^-2-0.5*nuhat^2*tauhat^-2)
                        *zdrawsd*thetadrawsd/tauhat));
  
  
  ##show(denominator)
  Z1mat<-matrix(rep(Z[,1],ntheta),ncol=ntheta);
  Z2mat<-matrix(rep(Z[,2],ntheta),ncol=ntheta);
  sigmamat<-matrix(rep(sigmaZ2,ntheta),ncol=ntheta);
  Thetadrawmat<-matrix(rep(Thetadraws[,1],n),nrow=n,byrow=TRUE)
  numerator<- apply(pztheta(Z1mat,Thetadrawmat)*exp(-0.5*(Z1mat-Thetadrawmat)^2
                                                    -0.5*tauhat^(-2)*(Thetadrawmat-nuhat)^(2) +0.5*(thetadrawsd^(-2))*Thetadrawmat^2)*thetadrawsd/tauhat * dnorm((Z2mat-Thetadrawmat)/sigmamat)/sigmamat,1,mean)
  
  ##show(numerator)
  logL=log(numerator/denominator);
 # %objective function; note the sign flip, since we are doing minimization
  
  LLH=-sum(logL); 
  return(list(LLH=LLH, logL=logL))
  
}

if (FALSE) {

nuhat<-5
tauhat<-1
betap<-c(0.3,0.4,0.5)
gammap<-c(0.3,0.4,0.5)
cutoffs<-c(0.6473313, 0.9394543)
symmetric<-1
set.seed(1)
Z<-matrix(rnorm(2*100),100,2)
Z<-cbind(Test$V1,Test$V2)
vcutoff<-1
set.seed(2)
Zdraws<-matrix(rnorm(2*50),50,2)
write.csv(Zdraws,"/Users/virasemenora/Dropbox (MIT)/Isaiah Andrews Replication/Zdraws.csv")

set.seed(3)
Thetadraws<-matrix(rnorm(2*25),25,2)
write.csv(Thetadraws,"/Users/virasemenora/Dropbox (MIT)/Isaiah Andrews Replication/Thetadraws.csv")

sigmaZ2<-rep(1,100)
zdrawsd <-1;
thetadrawsd<-1;

LLH<-ReplicationAnalyticLogLikelihood(nuhat, 
                                      tauhat, 
                                      betap, 
                                      cutoffs, 
                                      symmetric, 
                                      Z, 
                                      sigmaZ2)
}