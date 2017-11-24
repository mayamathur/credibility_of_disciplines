# not tested
source("ReplicationAnalyticLogLikelihood.R")
source("ReplicationpzthetaLogLikelihood.R")
source("VariationVarianceLogLikelihood.R")
source("SelectionTable.R")
source("RobustVariance.R")
#show(number)
MaxEval<-10^5
MaxIter<-10^5;
Tol<-10^(-8);
stepsize<-10^(-6);


if (identificationapproach==1) {
  
  
  if (symmetric == 0) {
    LLH <-function(Psi) {
      
      
      # ~~~~~ MM NOTE: THIS VERSION (SYMMETRIC == FALSE) ALLOWS MEAN TO BE NONZERO
      
      A<-ReplicationAnalyticLogLikelihood(Psi[1], Psi[2], c(Psi[-c(1,2)], 1),cutoffs,symmetric, Z, sigmaZ2);
      #A<-A$LLH
      return(A)
    }
  }
  
  
  else {
    if (spec_test == 1) {
      LLH <-function(Psi) {
       
        A<-ReplicationpzthetaLogLikelihood(0,Psi[1], c(Psi[2:(2+length(cutoffs)-1)], 1),
                                           c(rep (0, length(Psi[2:2+length(cutoffs)-1])),0),
                                           cutoffs,symmetric, Z, sigmaZ2,vcutoff, Zdraws, Thetadraws,zdrawsd,thetadrawsd);
        #A<-A$LLH
        return(A)
      }
      
      LLH_test <-function(Psi) {
        A<-ReplicationpzthetaLogLikelihood(0,Psi[1],  c(Psi[2:(2+length(cutoffs)-1)], 1),
                                           c(Psi[(2+length(cutoffs)):length(Psi)], 0),
                                           cutoffs,symmetric, Z, sigmaZ2,
                                           vcutoff, Zdraws, Thetadraws,zdrawsd,thetadrawsd);
        #A<-A$LLH
        return(A)
      } 
    }
    
    else {
      # ~~~~~ MM NOTE: THIS VERSION (SYMMETRIC == FALSE) FORCES THE MEAN TO BE 0
      
      LLH <-function (Psi) {
        A<-ReplicationAnalyticLogLikelihood(0,Psi[1], c(Psi[-1], 1),cutoffs,symmetric, Z, sigmaZ2);
        #A<-A$LLH
        return(A)
      }
    }
  }
  nn<<- number;
  
}




if (identificationapproach==2) {
  if (symmetric == 0) {
    
    if (symmetric_p == 1) {
      LLH <- function (Psi) {
        A<-VariationVarianceLogLikelihood(Psi[1], Psi[2], c(1, Psi[-c(1,2)], rev(c(Psi[-c(1,2)])), 1),cutoffs,symmetric, X[includeinestimation], sigma[includeinestimation],0);
        #A<-A$LLH
        return(A)
        }
      
      nn<<-sum(includeinestimation);
    } else {
    #  #show(c(Psi[-c(1,2)],  1))
      LLH <- function (Psi) {
        A<-VariationVarianceLogLikelihood(Psi[1], Psi[2], c(Psi[-c(1,2)],  1),cutoffs,symmetric, X[includeinestimation], sigma[includeinestimation],0);
        #A<-A$LLH
        return(A)
      }
      nn<<-sum(includeinestimation);
    }
    
  } else {
    LLH <-function (Psi) {
      A<-VariationVarianceLogLikelihood(0, Psi[1],c(Psi[-1], 1),cutoffs,symmetric, X[includeinestimation], sigma[includeinestimation],0);
      #A<-A$LLH
      return(A)
  
    }
    nn<<-sum(includeinestimation);  
  }
  #show(nn)
}


LLH_only<-function (Psi){
  A<-LLH(Psi);
  return(A$LLH)
}

if (application<=13) {
  #lower.b =rep(0,length(Psihat0))# why do we set lower range for thetabar to 0?
  lower.b = c(-Inf,rep(0,length(Psihat0)-1))
  
} else {
  lower.b = c(-Inf,rep(0,length(Psihat0)-1))
}

upper.b=rep(Inf,length(Psihat0))


#%find maximum likelihood estimator using just LLH
if (spec_test == 1) {
  

  findmin<-nlminb(objective=LLH_only, start=Psihat0,lower=lower.b,upper=upper.b);
  Psihat<-findmin$par;
  #show(Psihat)
  LLHmax<-findmin$objective;
  
  Psihat_test<-c(Psihat, rep(0,length(cutoffs)));
  score_mat<- matrix(0,number , length(Psihat_test));
  
  for (n1 in 1:length(Psihat_test)) {
    theta_plus<-Psihat_test;
    theta_plus[n1]<-theta_plus[n1]+stepsize;
    plus<-LLH_test(theta_plus);
    LLH_plus<-plus$LLH;
    logL_plus<-plus$logL;
 
    theta_minus<-Psihat_test;
    theta_minus[n1]<-theta_minus[n1]-stepsize;
    minus<-LLH_test(theta_minus);
    LLH_minus<-minus$LLH;
    logL_minus<-minus$logL;
    score_mat[,n1]<-(logL_plus-logL_minus)/(2*stepsize);
    
    
    
    
  }
  mean_score_mat <- as.matrix(apply(score_mat,2,mean))
  score_stat<-number*t(mean_score_mat)%*%solve(cov(score_mat))%*%mean_score_mat
  score_pval<-1 - pchisq(score_stat,length(cutoffs));

} else { # ~~~~ IF SPEC_TEST = 0 (THE WAY THEY SPECIFIED)

  findmin<-nlminb(objective=LLH_only, start=Psihat0,lower=lower.b,upper=upper.b,control = list(eval.max = MaxEval, iter.max = MaxIter, abs.tol = Tol));
 
  # ~~~~~~ MM NOTE: THIS PART IS IMPORTANT
  # I THINK FINDMIN IS THE ML ESTIMATES
  # INDEED, THE BELOW PSIHAT MATCHES ESTIMATES IN TABLE 2
  
  # IF I START NUHAT AT 0.2, I GET:
  # -0.006, 1.24, 0.008, 0.29
  # SO FAIRLY CLOSE TO WHAT THEY GOT IN TABLE 2
  # ALSO NEARLY 0 IF NUHAT IS SET TO 1
  
  Psihat<-findmin$par
  
  # ~~~~~~ MM ADDED
  cat("\n Psihat is: ", Psihat )
  
  
  LLHmax<-findmin$objective
  #Psihat<-findmin$estimate;
  #LLHmax<-findmin$minimum;
  #show(nn)
  if (application <=13) {
    Var_robust<-RobustVariance(stepsize, nn, Psihat, LLH,cluster_ID);
    se_robust<-sqrt(diag(Var_robust));
    
    SelectionTable(pathname,Psihat,se_robust,name,(identificationapproach==1),symmetric)
    
  }
 
}



# THIS PART GIVES ERROR MESSAGE FOR APPLICATION 3, BUT IT'S ONLY A SANITY CHECK, SO MAYBE OK?
# ONLY HAPPENS IF SYMMETRIC = 0
#% sanity check for replication studies
if (application <=13) {
  if (identificationapproach==1 && spec_test==0) {
    if (symmetric == 1) {
       LLH <- function (Psi) {
       A<-VariationVarianceLogLikelihood(0,Psi[1],c(Psi[-1], 1),cutoffs,symmetric,  X, sigma,1);
       #A<-A$LLH
       return(A)
      }
    }
    else {
      LLH <-function(Psi) {
        A<-VariationVarianceLogLikelihood(Psi[1], Psi[2], c(Psi[-c(1,2)], 1),cutoffs,symmetric, X, sigma,1);
        #A<-A$LLH
        return(A)
      }
    }
    LLH_only<-function (Psi){
      # #show(sigma)
      A<-LLH(Psi);
      return(A$LLH)
    }
    
    findmin<-nlminb(objective=LLH_only, start=Psihat0,lower=lower.b,upper=upper.b);
    PsihatSanity1<-findmin$par
    LLHmax<-findmin$objective
    #  PsihatSanity1<-findmin$estimate;
    #  LLHmax<-findmin$minimum;
    #  #show(PsihatSanity1)
    # #show(LLHmax)
    Var_robustSanity1<-RobustVariance(stepsize, nn, PsihatSanity1, LLH,cluster_ID);
    se_robustSanity1=sqrt(diag(Var_robustSanity1));
    #  #show(se_robustSanity1)
    SelectionTable(pathname,PsihatSanity1,se_robustSanity1, 'Sanitycheck',0,symmetric)
    
  } 
}
