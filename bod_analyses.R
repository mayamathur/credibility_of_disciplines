
############################## EXPLORE DATA ##############################



setwd("~/Dropbox/Personal computer/Independent studies/BOD (believability of disciplines)/bod_git")
s = read.csv("Data from Ioannidis/full_data.csv")
dim(s)  # expected: 23509

# 31% nonsignificant
# 6% marginally significant
# 63% significant
prop.table( table( s$t > 1.96 ) )
prop.table( table( s$t > 1.64 & s$t < 1.96 ) )
prop.table( table( s$t < 1.64 ) )


############################## RUN AK META-STUDY ALGORITHM ##############################

# run AK algorithm
setwd("~/Dropbox/Personal computer/Independent studies/BOD (believability of disciplines)/Analysis/AK code modified/Replication R Code/R Code")
# source("ApplicationScript.R")



############################## SIMULATE USING AK'S MLES ##############################

source("helper.R")

######## Z scale ########

# from running AK with Z-scores
bp0 = 1  # prob of publishing significant result (baseline)
bp1 = 0.3864854  
bp2 = 0.4464266
mu = 0
tau.tilde = 2.118972  

credibility( .n = 5000,
               .plot.n = 5000,
               .bp0 = bp0,
               .bp1 = bp1,
               .bp2 = bp2,
               .mu = mu,
               .tau.tilde = tau.tilde,
               .thresh = c(0, .1, .5),
               .incl.ref = TRUE,
               .scale = "Z",
               .plot.type = "ECDF.theta" )


######## X scale ########

# from running AK with raw effect sizes
bp0 = 1
bp1 = 0.441
bp2 = 0.491
tau = 0.584
mu = 0

# on X scale
credibility( .n = 5000,
             .plot.n = 5000,
             .bp0 = bp0,
             .bp1 = bp1,
             .bp2 = bp2,
             .mu = mu,
             .SEs = s$SE,
             .tau = tau,
             .thresh = c(0, .1, .5),
             .incl.ref = TRUE,
             .scale = "X",
             .plot.type = "ECDF.theta" )



