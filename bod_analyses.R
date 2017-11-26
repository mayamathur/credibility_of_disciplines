
############################## EXPLORE DATA ##############################

setwd("~/Dropbox/Personal computer/Independent studies/BOD (believability of disciplines)/bod_git")
s = read.csv("Data from Ioannidis/full_data.csv")
dim(s)  # expected: 23509



############################## RUN AK META-STUDY ALGORITHM ##############################

# run AK algorithm
setwd("AK code modified/Replication R Code/R")
# source("ApplicationScript.R")



############################## SIMULATE USING AK'S MLES ##############################


######## Z scale ########

# from running AK with Z-scores
bp0Z = 1  # prob of publishing significant result (baseline)
bp1Z = 0.3864854  
bp2Z = 0.4464266
mu = 0
tau.tilde = 2.118972  


setwd("~/Dropbox/Personal computer/Independent studies/BOD (believability of disciplines)/bod_git")

source("helper.R")
credibility( .n = 10000,
               .plot.n = 1500,
               .bp0 = bp0Z,
               .bp1 = bp1Z,
               .bp2 = bp2Z,
              # .bp0 = 1,
              # .bp1 = 0,
              # .bp2 = 0,
               .mu = mu,
               .tau.tilde = tau.tilde,
               .thresh = c(0, .5),
              .prox = 0.25,
               .incl.ref = TRUE,
               .scale = "Z",
               .plot.type = "scatter" )

# treating caliper as absolute deviation from true effect:
# without selectivity: 20%
# with MLE selectivity: 19%
# with horrible selectivity: 18%




######## X scale ########

# from running AK with raw effect sizes
bp0X = 1
bp1X = 0.441
bp2X = 0.491
tau = 0.584
mu = 0

setwd("~/Dropbox/Personal computer/Independent studies/BOD (believability of disciplines)/bod_git")
source("helper.R")
credibility( .n = 10000,
             .plot.n = 1000,
             # .bp0 = bp0X,
             # .bp1 = bp1X,
             # .bp2 = bp2X,
             .bp0 = 1,
             .bp1 = 0,
             .bp2 = 0,
             .mu = mu,
             .SEs = s$SE,
             .tau = tau,
             .thresh = c(0, 0.2),
             .prox = c(.1, .2),
             .incl.ref = TRUE,
             .scale = "X",
             .plot.type = "scatter" )



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                           RECONCILE WITH AK'S PLOT 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


dX = sim_X( .n = 100000, .mu = mu, .tau = tau, .SEs = s$SE * 4,
            .bp0 = bp0X, .bp1 = bp1X, .bp2 = bp2X ) 

# dX = sim_X( .n = 10000, .mu = mu, .tau = tau, .SEs = s$SE * 4,
#             .bp0 = 1, .bp1 = 0, .bp2 = 0 ) 



dZ = sim_Z( .n = 1000000, .mu = mu, .tau.tilde = tau.tilde,
            .bp0 = bp0Z, .bp1 = bp1Z, .bp2 = bp2Z ) 
dp = dZ[ dZ$publish == TRUE, ]


# distance from truth by significance
aggregate( abs( Xi - theta ) ~ abs( Zi ) > 1.96, data = dX, FUN = median)

# 0.005 threshold
aggregate( abs( Xi - theta ) ~ abs( Zi ) > qnorm( 1 - .005/2 ), data = dX, FUN = median)

# should match corrected inference (1.40 for median) in AK's generated CSV file
# run with n = 1,000,000
dp$close = abs( dp$Zi - 1.96 ) < .01
table(dp$close)


###### Sanity Check: Does AK's Corrected Inference Work? #####

# AK's correction
# look at the theta such that, conditional on that theta, the median observed Z is 1.96
# i.e., a horizontal cross-section of scatterplot

# check whether this holds up in the simulated data
# by considering the thetas that are close to the median-unbiased estimate
# each 
dp$close.theta = abs( dp$theta - 1.389 ) < .01
table(dp$close.theta)
prop.table( table( dp$Zi[dp$close.theta] < 1.96 ) )
# yes! close to 50% :)

dp$close.theta = abs( dp$theta - 1.594 ) < .01
table(dp$close.theta)
prop.table( table( dp$Zi[dp$close.theta] < 2.1 ) )
# yes! close to 50% :)

dp$close.theta = abs( dp$theta - 1.909 ) < .01
table(dp$close.theta)
prop.table( table( dp$Zi[dp$close.theta] < 2.31 ) )
# yes! close to 50% :)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                           RECONCILE WITH IOANNIDIS' ANALYSIS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# what are the true pre-study odds based on AK's MLEs and using Ioannidis' coarsening?
source("helper.R")

# 
# dZ = sim_Z( .n = 10000, .mu = mu, .tau.tilde = tau.tilde,
#            .bp0 = bp0Z, .bp1 = bp1Z, .bp2 = bp2Z ) 
# dZ = sim_Z( .n = 10000, .mu = mu, .tau.tilde = tau.tilde,
#             .bp0 = 1, .bp1 = 1, .bp2 = 1 ) 

dX = sim_X( .n = 10000, .mu = mu, .tau = tau, .SEs = s$SE,
                  .bp0 = bp0X, .bp1 = bp1X, .bp2 = bp2X ) 

# proportion true effects closer to d = 0 than d = 0.2
( tab = prop.table( table( abs( dX$Xi ) <= 0.1 ) ) )
( prob.H0 = tab[["TRUE"]] )

# correct coarsened prior odds of H0:H1
( prior.odds = prob.H0 / (1 - prob.H0) )
# only 14%!

# where on Ioannidis' X-axis is the truth?
log( prior.odds, base = 10 )

# try to reproduce his 13.5% 
frp_0( O = 25, alpha = 0.05, pwr = .5 )

frp_0( O = 1, alpha = 0.05, pwr = .5 )

















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




