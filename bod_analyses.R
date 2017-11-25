
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

<<<<<<< HEAD
=======
source("helper.R")

>>>>>>> 3327ed475b08aa40e5ed050cfabb1358f6b83fdb
######## Z scale ########

# from running AK with Z-scores
bp0Z = 1  # prob of publishing significant result (baseline)
bp1Z = 0.3864854  
bp2Z = 0.4464266
mu = 0
tau.tilde = 2.118972  

<<<<<<< HEAD
source("helper.R")
credibility( .n = 10000,
               .plot.n = 10000,
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

source("helper.R")
credibility( .n = 10000,
             .plot.n = 500,
             .bp0 = bp0X,
             .bp1 = bp1X,
             .bp2 = bp2X,
             # .bp0 = 1,
             # .bp1 = 0,
             # .bp2 = 0,
             .mu = mu,
             .SEs = s$SE,
             .tau = tau,
             .thresh = c( 0, 0.2 ),
             .prox = .1,
             .incl.ref = TRUE,
             .scale = "X",
             .plot.type = "dist" )

# median distance from truth:
# without selectivity:
# with MLE selectivity: 0.16
# with horrible selectivity: 0.17
# THIS SEEMS WRONG. FIX IT. 

# treating caliper as proportion of true effect: 
# without selectivity: 35% of published studies are within 25% of truth
# with MLE selectivity: 41%
# with horrible selectivity: 54%

# treating caliper as absolute deviation from true effect:
# without selectivity: 33%
# with MLE selectivity: 33%
# with horrible selectivity: 36%



# in scatterplot, being farther from 45-degree line => larger distance between est and truth


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                           RECONCILE WITH IOANNIDIS' ANALYSIS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# what are the true pre-study odds based on AK's MLEs and using Ioannidis' coarsening?
source("helper.R")

d = sim_X( .n = 10000, .mu = mu, .tau = tau, .SEs = s$SE,
                  .bp0 = bp0X, .bp1 = bp1X, .bp2 = bp2X ) 

# proportion true effects closer to d = 0 than d = 0.2
( tab = prop.table( table( abs( d$Xi ) <= 0.1 ) ) )
( prob.H0 = tab[["TRUE"]] )

# correct coarsened prior odds of H0:H1
( prior.odds = prob.H0 / (1 - prob.H0) )
# only 14%!

# where on Ioannidis' X-axis is the truth?
log( prior.odds, base = 10 )

# try to reproduce his 13.5% 
frp_0( O = 25, alpha = 0.05, pwr = .5 )

frp_0( O = 1, alpha = 0.05, pwr = .5 )

















=======
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

>>>>>>> 3327ed475b08aa40e5ed050cfabb1358f6b83fdb

######## X scale ########

# from running AK with raw effect sizes
bp0 = 1
bp1 = 0.441
bp2 = 0.491
tau = 0.584
mu = 0

<<<<<<< HEAD
=======
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
>>>>>>> 3327ed475b08aa40e5ed050cfabb1358f6b83fdb



