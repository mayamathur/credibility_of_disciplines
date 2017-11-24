library(ggplot2)
#library(grid)
library(gridExtra)
library(cowplot)

#rm(list=ls())
#application<-1;

##% whether we want to produce figures, and/or estimates, and/or bias corrections

Psihat.saved<-list()
##%simulations draws for p(z,theta) specifications
#pathname<-"/Users/virasemenora/Dropbox (MIT)/Isaiah Andrews Replication"
zdrawcount<-5*1000;
thetadrawcount<-5*1000;
zdrawsd<-5;
thetadrawsd<-4;
vcutoff<-1.96;
seed<-888
dofis1<-NULL
##%variable identificationapproach will be set later
##% 1: replications
##% 2: variation in variance

setup<-function(application) {
  show(application)
  if (application == 1) {

    
    filepath<-"/Applications/EconExperiments/sorted_names.csv";
    Studynames<-read.csv(paste0(pathname,filepath),header=FALSE)
    Studynames<<-Studynames$V1
    
    filepath<-"/Applications/EconExperiments/cleaned_econ_data.csv";
    my.data<-read.csv(paste0(pathname,filepath),header=FALSE)
    Z<-cbind(my.data[,1]/my.data[,2],my.data[,3]/my.data[,2])
    sigmaZ2<-my.data[,4]/my.data[,2];
    number=dim(Z)[1];
    cluster_ID=1:number;
    X=my.data[,1];
    sigma<-my.data[,2];
    identificationapproach=1;
    name="ReplicationEcon";
    #%Estimate baseline model, rather than running spec test
    spec_test=0;   
    #%Set number of normal draws for monte-carlo integration
    ndraws=10^4;
    #%Set cutoffs to use in step function: should be given in increasing order;
    cutoffs=1.96;
    #%Use a step function symmetric around zero
    symmetric=1;
    #%starting values for optimization
    #  #%[tau, betap(1)]
    Psihat0=c(2,1);
    
    args1<<-list(filepath=filepath,
                 my.data=my.data,
                 Z=Z,
                 sigmaZ2=sigmaZ2,
                 X=X,
                 number=number,
                 sigma=sigma,
                 cluster_ID=cluster_ID,
               Psihat0=Psihat0,
               name = name,
               symmetric= symmetric,
               cutoffs = cutoffs,
               spec_test = spec_test,
               identificationapproach = identificationapproach)
    #args1<<-c(3,3)
    
    
  } else if (application == 2) {
    filepath<-"/Applications/EconExperiments/sorted_names.csv";
    Studynames<-read.csv(paste0(pathname,filepath),header=FALSE)
    Studynames<<-Studynames$V1
    
    filepath<-"/Applications/EconExperiments/cleaned_econ_data.csv";
    my.data<-read.csv(paste0(pathname,filepath),header=FALSE)
    Z<-cbind(my.data[,1]/my.data[,2],my.data[,3]/my.data[,2])
    sigmaZ2<-my.data[,4]/my.data[,2];
    number=dim(Z)[1];
    cluster_ID=1:number;
    X=my.data[,1];
    sigma<-my.data[,2];
    
      
    
    X<-X/mean(sigma);
    sigma<-sigma/mean(sigma);
    identificationapproach=1;
    name="ReplicationEconSpecTest";
    #%Estimate baseline model, rather than running spec test
    spec_test=1;   
    #%Set number of normal draws for monte-carlo integration
    ndraws=10^4;
    #%Set cutoffs to use in step function: should be given in increasing order;
    cutoffs=1.96;
    #%Use a step function symmetric around zero
    symmetric=1;
    #%starting values for optimization
    #  #%[tau, betap(1)]
    Psihat0=c(1,1); 
    # 
    set.seed(seed)
    Zdraws<-as.matrix(rnorm(zdrawcount)*zdrawsd);
    Thetadraws<-as.matrix(rnorm(thetadrawcount)*thetadrawsd);
   # write.csv(Zdraws,paste0(pathname,"/Zdraws.csv"))
    #write.csv(Thetadraws,paste0(pathname,"/Thetadraws.csv"))
    Zdraws<-read.csv(paste0(pathname,"/Zdraws.csv"))
    Thetadraws<-read.csv(paste0(pathname,"/Thetadraws.csv"))
    Zdraws<-Zdraws[,2]
    Thetadraws<-Thetadraws[,2]
    args1<<-list(filepath=filepath,
                 Z=Z,
                 sigmaZ2=sigmaZ2,
                 X=X,
                 number=number,
                 Zdraws=Zdraws,
                 Thetadraws=Thetadraws,
                 sigma=sigma,
                 cluster_ID=cluster_ID,
                 Psihat0=Psihat0,
                 name = name,
                 symmetric= symmetric,
                 cutoffs = cutoffs,
                 spec_test = spec_test,
                 identificationapproach = identificationapproach,
                 my.data=my.data)
    
  }
  
  # ~~~~~~~ THE ONE WE'RE USING
  else if (application == 3) {
    
    Studynames<<-NULL
      
    filepath<-"/Applications/PsychExperiments/cleaned_psych_data.csv";
    my.data<-read.csv(paste0(pathname,filepath),header=FALSE);
    dofis1<-(my.data[,7]==1);
    use<-dofis1;
    my.data=my.data[use,]
    
    Z<-cbind(my.data[,1]/my.data[,2],my.data[,3]/my.data[,2])
    sigmaZ2<-my.data[,4]/my.data[,2];
    number=dim(Z)[1];
    cluster_ID=1:number;
    X=my.data[,1];
    sigma<-my.data[,2];
    identificationapproach=1;
    name="ReplicationPsych";
    #%Estimate baseline model, rather than running spec test
    spec_test=0;   
    #%Set number of normal draws for monte-carlo integration
    #ndraws=10^4;
    #%Set cutoffs to use in step function: should be given in increasing order;
    cutoffs=c(1.64,1.96);
    
    
    
    # ~~~~~~ MM NOTE: HERE IS WHERE THEY FORCE THE MEAN TO BE ZERO
    #%Use a step function symmetric around zero
    symmetric=1;
    #%starting values for optimization
    #  #%[tau, betap(1)]
    Psihat0=c(1,1,1);
    
    # ~~~ MM CHANGED USING BELOW PART TO ALLOW NONZERO MEAN
    # SEE ESTIMATINGSELECTION.R
    # I BELIEVE THESE ARE JUST THE STARTING VALUES FOR OPTIMIZATION PER THE ABOVE,
    # SO LET'S JUST START AT MEAN 0
#     symmetric=0
#     Psihat0 = c(0,1,1,1)
    
    
    args1<<-list(filepath=filepath,
                 dofis1=dofis1,
                 Z=Z,
                 sigmaZ2=sigmaZ2,
                 X=X,
                 number=number,
                 sigma=sigma,
                 cluster_ID=cluster_ID,
                 Psihat0=Psihat0,
                 name = name,
                 symmetric= symmetric,
                 cutoffs = cutoffs,
                 spec_test = spec_test,
                 identificationapproach = identificationapproach,
                 my.data=my.data)
    
    
  } else if (application == 4) {
    Studynames<<-NULL
    
   # filepath<-"/Applications/PsychExperiments/sorted_1stAuthor.csv";
   # Studynames<-read.csv(paste0(pathname,filepath),header=TRUE)
   # Studynames<-Studynames$x1st_author
    filepath<-"/Applications/PsychExperiments/cleaned_psych_data.csv";
    my.data<-read.csv(paste0(pathname,filepath),header=FALSE);
    dofis1<-(my.data[,7]==1);
    approved <-as.logical(my.data[,5])
    use<-dofis1 & approved;
    my.data=my.data[use,]
    
    Z<-cbind(my.data[,1]/my.data[,2],my.data[,3]/my.data[,2])
    sigmaZ2<-my.data[,4]/my.data[,2];
    number=dim(Z)[1];
    cluster_ID=1:number;
    X=my.data[,1];
    sigma<-my.data[,2];
    identificationapproach=1;
    name="ReplicationPsychApproved";
    #%Estimate baseline model, rather than running spec test
    spec_test=0;   
    #%Set number of normal draws for monte-carlo integration
    #ndraws=10^4;
    #%Set cutoffs to use in step function: should be given in increasing order;
    cutoffs=c(1.64,1.96);
    #%Use a step function symmetric around zero
    symmetric=1;
    #%starting values for optimization
    #  #%[tau, betap(1)]
    Psihat0=c(1,1,1);
    
    args1<<-list(filepath=filepath,
                 dofis1=dofis1,
                 Z=Z,
                 sigmaZ2=sigmaZ2,
                 X=X,
                 number=number,
                 sigma=sigma,
                 cluster_ID=cluster_ID,
                 Psihat0=Psihat0,
                 name = name,
                 symmetric= symmetric,
                 cutoffs = cutoffs,
                 spec_test = spec_test,
                 identificationapproach = identificationapproach,
                 my.data=my.data)
  } else if (application == 5) {
    filepath<-"/Applications/PsychExperiments/sorted_1stAuthor.csv";
    Studynames<-read.csv(paste0(pathname,filepath),header=TRUE)
    NAMES<<-as.matrix(Studynames$x1st_author)
    
    filepath<-"/Applications/PsychExperiments/cleaned_psych_data.csv";
    my.data<-read.csv(paste0(pathname,filepath),header=FALSE);
    dofis1<-(my.data[,7]==1);
    #approved <-as.logical(my.data[,5])
    denominatordof30<-my.data[,8]==0|my.data[,8]>=30;
    use<-dofis1 &denominatordof30;
    my.data=my.data[use,]
    show(use)
    show(Studynames)
    show(class(Studynames))
    NAMES<<-NAMES[use]
    Studynames<<-NAMES
    Studynames<-NULL
    
    Z<-cbind(my.data[,1]/my.data[,2],my.data[,3]/my.data[,2])
    sigmaZ2<-my.data[,4]/my.data[,2];
    number=dim(Z)[1];
    cluster_ID=1:number;
    X=my.data[,1];
    sigma<-my.data[,2];
    identificationapproach=1;
    name="ReplicationPsychLargedof";
    #%Estimate baseline model, rather than running spec test
    spec_test=0;   
    #%Set number of normal draws for monte-carlo integration
    #ndraws=10^4;
    #%Set cutoffs to use in step function: should be given in increasing order;
    cutoffs=c(1.64,1.96);
    #%Use a step function symmetric around zero
    symmetric=1;
    #%starting values for optimization
    #  #%[tau, betap(1)]
    Psihat0=c(2,1,1);
    
    args1<<-list(filepath=filepath,
                 dofis1=dofis1,
                 Z=Z,
                 sigmaZ2=sigmaZ2,
                 X=X,
                 number=number,
                 sigma=sigma,
                 cluster_ID=cluster_ID,
                 Psihat0=Psihat0,
                 name = name,
                 symmetric= symmetric,
                 cutoffs = cutoffs,
                 spec_test = spec_test,
                 identificationapproach = identificationapproach,
                 my.data=my.data)
    
  } else if (application == 6) {
    
    filepath<-"/Applications/PsychExperiments/sorted_1stAuthor.csv";
    Studynames<-read.csv(paste0(pathname,filepath),header=TRUE)
    Studynames<<-Studynames$x1st_author
    Studynames<-NULL
    
    filepath<-"/Applications/PsychExperiments/cleaned_psych_data.csv";
    my.data<-read.csv(paste0(pathname,filepath),header=FALSE)
    Z<-cbind(my.data[,1]/my.data[,2],my.data[,3]/my.data[,2])
    sigmaZ2<-my.data[,4]/my.data[,2];
    number=dim(Z)[1];
    cluster_ID=1:number;
    X=my.data[,1];
    sigma<-my.data[,2];
    
    
    
    X<-X/mean(sigma);
    sigma<-sigma/mean(sigma);
    identificationapproach=1;
    name="ReplicationPsychSpecTest";
    #%Estimate baseline model, rather than running spec test
    spec_test=1;   
    #%Set number of normal draws for monte-carlo integration
    ndraws=10^4;
    #%Set cutoffs to use in step function: should be given in increasing order;
    cutoffs=c(1.64,1.96);
    #%Use a step function symmetric around zero
    symmetric=1;
    #%starting values for optimization
    #  #%[tau, betap(1)]
    Psihat0=c(1.4698,0.02,0.2843); 
    # 
    set.seed(seed)
    Zdraws<-as.matrix(rnorm(zdrawcount)*zdrawsd);
    Thetadraws<-as.matrix(rnorm(thetadrawcount)*thetadrawsd);
    
    Zdraws<-read.csv(paste0(pathname,"/Zdraws.csv"))
    Thetadraws<-read.csv(paste0(pathname,"/Thetadraws.csv"))
    Zdraws<-Zdraws[,2]
    Thetadraws<-Thetadraws[,2]
    args1<<-list(filepath=filepath,
                 Z=Z,
                 sigmaZ2=sigmaZ2,
                 X=X,
                 number=number,
                 Zdraws=Zdraws,
                 Thetadraws=Thetadraws,
                 sigma=sigma,
                 cluster_ID=cluster_ID,
                 Psihat0=Psihat0,
                 name = name,
                 symmetric= symmetric,
                 cutoffs = cutoffs,
                 spec_test = spec_test,
                 identificationapproach = identificationapproach,
                 my.data=my.data)
    
  } else if (application == 7) {
    
    filepath<-"/Applications/MinimumWagev2/cleaned_minwage_data.csv";
    my.data<-read.csv(paste0(pathname,filepath),header=FALSE);
   
    Studynames<<-NULL
    
    X=my.data[,1];
    X=-X;
    number=length(X);
    cluster_ID=my.data[,3];
    
    includeinfigure<-as.logical(rep(1,number));
    includeinestimation<-as.logical(rep(1,number))
    
    sigma<-my.data[,2];
    identificationapproach=2;
    name="MinimumWage";
    #%Estimate baseline model, rather than running spec test
    spec_test=0;   
    #%Set number of normal draws for monte-carlo integration
    #ndraws=10^4;
    #%Set cutoffs to use in step function: should be given in increasing order;
    cutoffs=c(-1.96,0,1.96);
    #%Use a step function symmetric around zero
    symmetric=0;
    symmetric_p=0;
    #%starting values for optimization
    #  #%[tau, betap(1)]
    Psihat0=c(0,1.1,1,1,1);
    
    args1<<-list(filepath=filepath,
                 dofis1=dofis1,
                 Z=Z,
                 sigmaZ2=sigmaZ2,
                 X=X,
                 number=number,
                 sigma=sigma,
                 cluster_ID=cluster_ID,
                 Psihat0=Psihat0,
                 name = name,
                 symmetric= symmetric,
                 symmetric_p=symmetric_p,
                 cutoffs = cutoffs,
                 spec_test = spec_test,
                 identificationapproach = identificationapproach,
                 includeinestimation=includeinestimation,
                 my.data=my.data)
  }else if (application == 8) {
    filepath<-"/Applications/MinimumWagev2/cleaned_minwage_data.csv";
    my.data<-read.csv(paste0(pathname,filepath),header=FALSE);
    
    Studynames<<-NULL
    
    X=my.data[,1];
    X=-X;
    sigma<-my.data[,2];
    cluster_ID=my.data[,3];
    
    published<-my.data[,4];
    use = published==1;
    X = X[use];
    sigma = sigma[use];
    cluster_ID=cluster_ID[use];
    
    number=length(X);
    includeinfigure<-as.logical(rep(1,number));
    includeinestimation<-as.logical(rep(1,number))
    
    
    identificationapproach=2;
    name="MinimumWagePublished";
    #%Estimate baseline model, rather than running spec test
    spec_test=0;   
    #%Set number of normal draws for monte-carlo integration
    #ndraws=10^4;
    #%Set cutoffs to use in step function: should be given in increasing order;
    cutoffs=c(-1.96,0,1.96);
    #%Use a step function symmetric around zero
    symmetric=0;
    symmetric_p=0;
    #%starting values for optimization
    #  #%[tau, betap(1)]
    Psihat0=c(0,1,1,1,1);
    
    args1<<-list(filepath=filepath,
                 dofis1=dofis1,
                 Z=Z,
                 sigmaZ2=sigmaZ2,
                 X=X,
                 number=number,
                 sigma=sigma,
                 cluster_ID=cluster_ID,
                 Psihat0=Psihat0,
                 name = name,
                 symmetric= symmetric,
                 symmetric_p=symmetric_p,
                 cutoffs = cutoffs,
                 spec_test = spec_test,
                 identificationapproach = identificationapproach,
                 includeinestimation=includeinestimation,
                 my.data=my.data)
    
  } else if (application == 9) {
    filepath<-"/Applications/MinimumWagev2/cleaned_minwage_data.csv";
    my.data<-read.csv(paste0(pathname,filepath),header=FALSE);
    
    Studynames<<-NULL
    
    X=my.data[,1];
    X=-X;
    sigma<-my.data[,2];
    cluster_ID=my.data[,3];
    
    published<-my.data[,4];
    use = published==1;
    X = X[use];
    sigma = sigma[use];
    cluster_ID=cluster_ID[use];
    
    number=length(X);
    includeinfigure<-as.logical(rep(1,number));
    includeinestimation<-as.logical(rep(1,number))
    
    
    identificationapproach=2;
    name="MinimumWageNoSelection";
    #%Estimate baseline model, rather than running spec test
    spec_test=0;   
    #%Set number of normal draws for monte-carlo integration
    #ndraws=10^4;
    #%Set cutoffs to use in step function: should be given in increasing order;
    cutoffs=c(0);
    #%Use a step function symmetric around zero
    symmetric=0;
    symmetric_p=1;
    #%starting values for optimization
    #  #%[tau, betap(1)]
    Psihat0=c(0,1);
    
    args1<<-list(filepath=filepath,
                 dofis1=dofis1,
                 Z=Z,
                 sigmaZ2=sigmaZ2,
                 X=X,
                 number=number,
                 sigma=sigma,
                 cluster_ID=cluster_ID,
                 Psihat0=Psihat0,
                 name = name,
                 symmetric= symmetric,
                 symmetric_p=symmetric_p,
                 cutoffs = cutoffs,
                 spec_test = spec_test,
                 identificationapproach = identificationapproach,
                 includeinestimation=includeinestimation,
                 my.data=my.data)
    
  }else if (application == 10) {
    filepath<-"/Applications/deworming/cleaned_deworming_data.csv";
    my.data<-read.csv(paste0(pathname,filepath),header=FALSE);
    
    filepath<-"/Applications/deworming/DewormingLabels.csv";
    Studynames<-read.csv(paste0(pathname,filepath),header=TRUE)
    Studynames<<-Studynames$`Shorted.Name..year`
    Studynames<-NULL
    
    
    X=my.data[,1];
    sigma<-my.data[,2];
    cluster_ID=my.data[,3];
    number=length(X);
    
   
    
   
    includeinfigure<-as.logical(rep(1,number));
    includeinestimation<-as.logical(rep(1,number))
    
    
    identificationapproach=2;
    name="Deworming";
    #%Estimate baseline model, rather than running spec test
    spec_test=0;   
    #%Set number of normal draws for monte-carlo integration
    #ndraws=10^4;
    #%Set cutoffs to use in step function: should be given in increasing order;
    cutoffs=c(-1.96,0,1.96);
    #%Use a step function symmetric around zero
    symmetric=0;
    symmetric_p=1;
    #%starting values for optimization
    #  #%[tau, betap(1)]
    Psihat0=c(0,1,1);
    
    args1<<-list(filepath=filepath,
                 
                 dofis1=dofis1,
                 Z=Z,
                 sigmaZ2=sigmaZ2,
                 X=X,
                 number=number,
                 sigma=sigma,
                 cluster_ID=cluster_ID,
                 Psihat0=Psihat0,
                 name = name,
                 symmetric= symmetric,
                 symmetric_p = symmetric_p,
                 cutoffs = cutoffs,
                 spec_test = spec_test,
                 identificationapproach = identificationapproach,
                 includeinestimation=includeinestimation,
                 my.data=my.data)
  } else if (application == 11) {
    filepath<-"/Applications/deworming/cleaned_deworming_data.csv";
    my.data<-read.csv(paste0(pathname,filepath),header=FALSE);
    
    filepath<-"/Applications/deworming/DewormingLabels.csv";
    Studynames<-read.csv(paste0(pathname,filepath),header=TRUE)
    Studynames<<-Studynames$`Shorted.Name..year`
    Studynames<-NULL
    X=my.data[,1];
    sigma<-my.data[,2];
    cluster_ID=my.data[,3];
    number=length(X);
    
    
    
    
    includeinfigure<-as.logical(rep(1,number));
    includeinestimation<-as.logical(rep(1,number))
    
    
    identificationapproach=2;
    name="DewormingAsymmetric";
    #%Estimate baseline model, rather than running spec test
    spec_test=0;   
    #%Set number of normal draws for monte-carlo integration
    #ndraws=10^4;
    #%Set cutoffs to use in step function: should be given in increasing order;
    cutoffs=c(-1.96,0,1.96);
    #%Use a step function symmetric around zero
    symmetric=0;
    symmetric_p=0;
    #%starting values for optimization
    #  #%[tau, betap(1)]
    Psihat0=c(0,1,1,1,1);
    
    args1<<-list(filepath=filepath,
                 dofis1=dofis1,
                 Z=Z,
                 sigmaZ2=sigmaZ2,
                 X=X,
                 number=number,
                 sigma=sigma,
                 cluster_ID=cluster_ID,
                 Psihat0=Psihat0,
                 name = name,
                 symmetric= symmetric,
                 symmetric_p = symmetric_p,
                 cutoffs = cutoffs,
                 spec_test = spec_test,
                 identificationapproach = identificationapproach,
                 includeinestimation=includeinestimation,
                 my.data=my.data)
    
  } else if (application == 12) {
    filepath<-"/Applications/deworming/cleaned_deworming_data.csv";
    my.data<-read.csv(paste0(pathname,filepath),header=FALSE);
    
    filepath<-"/Applications/deworming/DewormingLabels.csv";
    Studynames<-read.csv(paste0(pathname,filepath),header=TRUE)
    Studynames<<-Studynames$`Shorted.Name..year`
    Studynames<-NULL
    X=my.data[,1];
    sigma<-my.data[,2];
    cluster_ID=my.data[,3];
    number=length(X);
    
    
    
    
    includeinfigure<-as.logical(rep(1,number));
    includeinestimation<-as.logical(rep(1,number))
    
    
    identificationapproach=2;
    name="DewormingRestrictedAsymmetric";
    #%Estimate baseline model, rather than running spec test
    spec_test=0;   
    #%Set number of normal draws for monte-carlo integration
    #ndraws=10^4;
    #%Set cutoffs to use in step function: should be given in increasing order;
    cutoffs=c(0);
    #%Use a step function symmetric around zero
    symmetric=0;
    symmetric_p=0;
    #%starting values for optimization
    #  #%[tau, betap(1)]
    Psihat0=c(0,1,1);
    
    args1<<-list(filepath=filepath,
                 dofis1=dofis1,
                 Z=Z,
                 sigmaZ2=sigmaZ2,
                 X=X,
                 number=number,
                 sigma=sigma,
                 cluster_ID=cluster_ID,
                 Psihat0=Psihat0,
                 name = name,
                 symmetric= symmetric,
                 symmetric_p = symmetric_p,
                 cutoffs = cutoffs,
                 spec_test = spec_test,
                 identificationapproach = identificationapproach,
                 includeinestimation=includeinestimation,
                 my.data=my.data)
  } else if (application == 13) {
    filepath<-"/Applications/deworming/cleaned_deworming_data.csv";
    my.data<-read.csv(paste0(pathname,filepath),header=FALSE);
    
    filepath<-"/Applications/deworming/DewormingLabels.csv";
    Studynames<-read.csv(paste0(pathname,filepath),header=TRUE)
    Studynames<<-Studynames$`Shorted.Name..year`
    Studynames<-NULL
    X=my.data[,1];
    sigma<-my.data[,2];
    cluster_ID=my.data[,3];
    number=length(X);
    
    
    
    
    includeinfigure<-as.logical(rep(1,number));
    includeinestimation<-as.logical(rep(1,number))
    
    
    identificationapproach=2;
    name="DewormingNoSelection";
    #%Estimate baseline model, rather than running spec test
    spec_test=0;   
    #%Set number of normal draws for monte-carlo integration
    #ndraws=10^4;
    #%Set cutoffs to use in step function: should be given in increasing order;
    cutoffs=c(0);
    #%Use a step function symmetric around zero
    symmetric=0;
    symmetric_p=1;
    #%starting values for optimization
    #  #%[tau, betap(1)]
    Psihat0=c(0,1);
    
    args1<<-list(filepath=filepath,
                 dofis1=dofis1,
                 Z=Z,
                 sigmaZ2=sigmaZ2,
                 X=X,
                 number=number,
                 sigma=sigma,
                 cluster_ID=cluster_ID,
                 Psihat0=Psihat0,
                 name = name,
                 symmetric= symmetric,
                 symmetric_p = symmetric_p,
                 cutoffs = cutoffs,
                 spec_test = spec_test,
                 identificationapproach = identificationapproach,
                 includeinestimation=includeinestimation,
                 my.data=my.data)
  }
  else if (application == 14) {
    my.data<-NULL
    set.seed(888)
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
    
 
    number=dim(Z)[1];
    cluster_ID=1:number;
    X=my.data[,1];
    sigma<-my.data[,2];
    identificationapproach=1;
    name="Monte Carlo";
    #%Estimate baseline model, rather than running spec test
    spec_test=0;   
    #%Set number of normal draws for monte-carlo integration
    ndraws=10^4;
    #%Set cutoffs to use in step function: should be given in increasing order;
    cutoffs=1.96;
    #%Use a step function symmetric around zero
    symmetric=0;
    #%starting values for optimization
    #  #%[tau, betap(1)]
    Psihat0=c(0,2,1);
    
    args1<<-list(filepath=filepath,
                 my.data=my.data,
                 Z=Z,
                 sigmaZ2=sigmaZ2,
                 X=X,
                 number=number,
                 sigma=sigma,
                 cluster_ID=cluster_ID,
                 Psihat0=Psihat0,
                 name = name,
                 symmetric= symmetric,
                 cutoffs = cutoffs,
                 spec_test = spec_test,
                 identificationapproach = identificationapproach)
  }
  
  
  return(args1)
}



