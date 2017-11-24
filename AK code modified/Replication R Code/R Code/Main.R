# clear memory
rm(list=ls())
# insert your path name here
pathname = "~/Dropbox/Personal computer/Independent studies/BOD (believability of disciplines)/Code from Andrews and Kasy/Replication R Code"
#pathname<-"/Users/virasemenora/Dropbox (MIT)/Isaiah Andrews Replication"
#pathname<-"~/Dropbox (MIT)/2016Summer/PublicationBias/Replication R Code"
setwd(pathname)
# Create a folder "FiguresandTables" under this path to save figures and tables
setwd("R Code")

#install.packages("ggplot2","grid","gridExtra","cowplot","reshape")
#install.packages(c("ggplot2","gridExtra","reshape","mvtnorm","cowplot")) #grid is a base package now, cowplot needs >= 3.3.0
source("AllApplications.R")
# set numbers of applications to run

###%##%
#  ##%picking which application to run
###% 1: econ experiments
###% 2: econ experiments, p(z,theta) specification test
###% 3: psychology experiments
###% 4: psychology experiments, approved replications
###% 5: psychology experiments, denominator dof><-30
###% 6: psychology experiments, p(z,theta) specification test
###% 7: minimum wage effect
###% 8: minimum wage effect, published
###% 9: minimum wage, no selection
###% 10: deworming
###% 11: deworming, asymmetric
###% 12: deworming, asymmetric, single discontinuity
###% 13: deworming, no selection

#applications<-c(1,3,4,5,7,8,9,10,11,12,13)
applications<-c(3)
# Choose operations to do
dofigures = 1;
doestimates = 1;
docorrections = 1;

# MM note: application 3 creates Table 2 

# MM, to allow stepping through below code
application = applications

for (application in applications) {
  application<<-application
  show(application)
  #  show(application)
  
  args1<<-setup(application)
  attach(args1)
  
  if (doestimates) {
    source("EstimatingSelection.R")
    Psihat.saved[[application]]<-Psihat
    show(application)
    show(Psihat)
  }
  #show(application)
  if (dofigures) {
    critval<<-1.96
    shape<<-21
    source("DescriptiveStats.R")
    source("RDescriptiveStats2.R")
  
    if (identificationapproach==1) {
       
      source("RDescriptiveStats1.R")
      #  xaxisname<<-"Z"
      #  yaxisname<<-"ZRep"
      RDescriptiveStats1(Z,NULL,NULL,FALSE,name,symmetric)
      RDescriptiveStats2(my.data[,1],my.data[,2],paste0(name,'Sanitycheck'),symmetric,cluster_ID)
      # xaxisname<<-"W"
      #  yaxisname<<-"WRep"
      RDescriptiveStats1(Z,my.data[,1],my.data[,2],TRUE,name,symmetric)
      
      
    } else {
      
      RDescriptiveStats2(X,sigma,name,symmetric,cluster_ID)
      
    }
  }
  if (docorrections) {
    source("HorizontalBars.R")
    
    
  }
  
}
