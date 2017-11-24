#this is a wrapper script to estimate the selection model and calculate
#selection corrected estimates based on Andrews and Kasy (2017)

# clear memory
#rm(list=ls())
# insert your path name here
pathname = "~/Dropbox/Personal computer/Independent studies/BOD (believability of disciplines)/Analysis/AK code modified/Replication R Code"
setwd(pathname)
# Create a folder "FiguresandTables" under this path to save figures and tables
setwd("R Code")

#install.packages(c("ggplot2","gridExtra","reshape","mvtnorm","cowplot")) # cowplot needs >= 3.3.0
library(ggplot2)
library(gridExtra)
library(cowplot)

# To replicate results in the main text, refer to AllApplicationsm and set the variable application accordingly. Set it to zero unless you are replicating results in the main text
#  (and in that case, you should be using AllApplications.m
application <- 0;
# Estimate baseline model so set it to 0, rather than running spec test against p(z, theta)
spec_test=0;   

##### 
#data preparation and parametric specification
#Based on your study type (replication studies or meta studies), edit and
#run corresponding section below for data preparation.
#identificationapproach = 1 if replication studies;
#identificationapproach = 2 if meta studies.

identificationapproach=2;


if (identificationapproach == 2) {
      # data is formatted as [original study X1, sigma1] 
  
      # random sample of data
      # filepath = "~/Dropbox/Personal computer/Independent studies/BOD (believability of disciplines)/Analysis/Data from Ioannidis/random_sample.csv"
      # 
      # full data
      filepath = "~/Dropbox/Personal computer/Independent studies/BOD (believability of disciplines)/Analysis/Data from Ioannidis/full_data.csv"
      
      my.data<-read.csv(filepath,header=TRUE)
      
      # exclude 1 study with no SE
      my.data = my.data[ !is.na(my.data$SE), ]
      #my.data<-read.csv(filepath,header=FALSE);
      
      Studynames<<-NULL
      
      #Name for outputs
      name="Ioannidis";
      
      #Assign these variables based on your dataset and parametric specification but do not change the variable names
      # MM: treating t-values as z-scores
      
      # TO DO ANALYSIS ON EFFECT-SIZE SCALE (TAU)
      X = my.data$D
      sigma = my.data$SE
      
      # TO DO ANALYSIS ON Z-SCORE SCALE (TAU-TILDE)
      # X = my.data$tr
      # sigma = rep( 1, length(X) )
      
      # from the paper, this is an ID for study
      #  because for min wage data, some studies contributed multiple point estimates
      #  I will treat as independent
      cluster_ID = 1:nrow(my.data)
      
      #indicator variable 1=the estimate is used in estimation of the selection model - see EstimatingSelection.R
      includeinestimation<-as.logical(rep(1,length(X)));

      # Set cutoffs to use in step function: should be given in increasing order;
      
      #cutoffs=c(1.64,1.96) # 2-cutoff case
      cutoffs = c(1.64, 1.96) # 5-cutoff case
      
      # Use a step function symmetric around zero
      symmetric=1
     # symmetric_p=0;
      # starting values for optimization
      # [thetabar, tau, betap(i)];dimension of betap(i) corresponds to cutoff
      #Psihat0=c( 0, 1.1, 1, 1, 1 );
      
      # TAKEN FROM PSYCH SCRIPT
      Psihat0=c( 1, 1, 1 ) # 2-cutoff case
      #Psihat0 = c( 1, 1, 1, 1, 1, 1 ) # 5-cutoff case
      
      number=length(X);
}



#####
#producing figures: 
#binned density plot for Z1; 
#scatter plot of initial Z1 against replication Z2 (for replication studies only); 
#scatter plot of initial X1 against standard errors (for replication studies and meta studies)
  critval<<-1.96 #set the cutoff displayed in figures of study estimates
  shape<<-21 #set the shape in ggplot scatter plots of study estimates 21=filled circle

  source("RDescriptiveStats2.R")
  
  # if (identificationapproach==1) {
  #   
  #   source("RDescriptiveStats1.R")
  # 
  #   RDescriptiveStats1(Z,NULL,NULL,FALSE,name,symmetric)
  # 
  #   RDescriptiveStats2(X,sigma,paste0(name,'Sanitycheck'),symmetric,cluster_ID)
  #   
  #   RDescriptiveStats1(Z,X,sigma,TRUE,name,symmetric)
  #   
  #   
  # } else {
  #   
  #   RDescriptiveStats2(X,sigma,name,symmetric,cluster_ID)
  #   
  # }

  
#####
#estimating the model
  
  source("EstimatingSelection.R")

  show(Psihat)


#####
#producing bias-corrected estimates and confidence sets
#  source("HorizontalBars.R")
  
#output selection bias corrected estimates and confidence sets to a csv file:
#Z1 is the original z-statistic, Z1_M is median unbiased estimator, Z1_L and Z1_U
#form a selection corrected confidence set for the estimate.
AdjustedZtable = cbind(as.character(Studynames), Z1, Z1_L, Z1_M, Z1_U )
write.csv(AdjustedZtable, file =paste0('../FiguresandTables/', name, 'AdjustedZtable.csv'),row.names = F)
  

