
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
               # .bp0 = bp0Z,
               # .bp1 = bp1Z,
               # .bp2 = bp2Z,
              .bp0 = 1,
              .bp1 = 0,
              .bp2 = 0,
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
             .plot.n = 10000,
             # .bp0 = bp0X,
             # .bp1 = bp1X,
             # .bp2 = bp2X,
             .bp0 = 1,
             .bp1 = 0,
             .bp2 = 0,
             .mu = mu,
             .SEs = s$SE*3,
             .tau = tau,
             .thresh = c(0, 0.2),
             .prox = c(.1, .5),
             .incl.ref = TRUE,
             .scale = "X",
             .plot.type = "scatter" )

# so with selectivity, Z > 1.96 findings aren't any farther from the truth than without selectivity
# but they tend to be biased by being LARGER, not smaller,
# so there is newfound bias. 
# this explains the proportion larger magnitude than truth as well. 



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                           RECONCILE WITH JM's ANALYSIS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# "On the Reproducibility of Psychological Science"

######## Empirical Data ########
# get RPP data for comparison
setwd("~/Dropbox/Personal computer/Miscellaneous/My publications/Published papers/OSC Estimating the reproducibility of psychological science/rpp_reproduce")
r = read.csv( "rpp_final_data.csv", header = TRUE )
# sanity check: should be 73 studies
nrow(r)

# randomly assign signs
set.seed(411); sign = ifelse( rbinom( n = nrow(r), prob = 0.5, size = 1 ) == 1, 1, -1 )
r$Zi = ( r$fis.o / r$sei.o ) * sign
table( r$Zi > 0 )



######## Simulate Under AK's Model ########
setwd("~/Dropbox/Personal computer/Independent studies/BOD (believability of disciplines)/bod_git")
source("helper.R")

AK = sim_X( .n = 10000, .mu = 0, .tau = 0.252, .SEs = r$sei.o,
            .bp0 = 1, .bp1 = 0.025, .bp2 = 0.375 ) 
AKp = AK[ AK$publish == TRUE, ]

prop.table( table( abs(AKp$theta) < 0.1 ) )


######## Simulate Under JM's Normal Model ########
# simulate true effect sizes and observed ES under their fitted normal model
#   which they say does not fit very well

JM1 = sim_X_JM( .n = 10000, .SEs = r$sei.o, .model = "normal" )
JM1p = JM1[ JM1$publish == TRUE, ]

# look at the point mass
#hist( JM1$theta, breaks = 20 )
# whoa


prop.table( table( abs(JM1p$theta) < 0.1 ) )



######## Simulate Under JM's Bimodal Moment Model ########
# simulate true effect sizes and observed ES under their fitted normal model
#   which they say does not fit very well

JM2 = sim_X_JM( .n = 10000, .mu = 0, .tau = 0.0877, .pi0 = 0.930, .SEs = r$sei.o,
                .alpha = 0.00569, .model = "moment" )
JM2p = JM2[ JM2$publish == TRUE, ]


# look at the point mass
hist( JM2$theta, breaks = 20 )
# whoa

# what proportion are < d = 0.1? (i.e., Ioannidis' coarsening)
prop.table( table( abs(JM2p$theta) < 0.1 ) )

# BOOKMARK
( prior.odds = .49 / (1 - .49) )


######## Competition! ########

colors = c("blue", "black", "red", "orange")

E2 = ggplot( data = AKp, aes( AKp$Zi ) ) +
  
  # stat_ecdf( data = AK, aes( x = AK$Zi, color = "No selectivity (publish everything)" ),
  #            show.legend = TRUE, lwd = 1 ) +  # ECDF without selectivity (all findings)
  
  stat_ecdf( aes( color = "AK model" ),
             show.legend = TRUE, lwd = 1.5 ) +  # ECDF of simulated data
  
  stat_ecdf( data = JM1p, aes( x = JM1p$Zi, color = "JM (normal model)" ),
             show.legend = TRUE, lwd = 1 ) + 
  
  stat_ecdf( data = JM2p, aes( x = JM2p$Zi, color = "JM (moment model)" ),
             show.legend = TRUE, lwd = 1 ) + 
  
  stat_ecdf( data = r, aes( x = r$Zi, color = "Empirical (random signs)" ),
             show.legend = TRUE, lwd = 1 ) +  # ECDF without selectivity (all findings)
  
  theme_bw() +
  scale_color_manual( name = " ", values = colors) +
  scale_x_continuous(limits=c(-6,6), breaks=seq(-6, 6, 1)) +
  scale_y_continuous(limits=c(0,1), breaks=seq(0, 1, .05)) +
  geom_vline( xintercept = c(-1.64, 1.64, 1.96, -1.96), linetype=2, color = "red" ) +
  xlab("Z-score") + ylab("Proportion below") +
  ggtitle("Model fit: ECDF of published Z-scores")
plot(E2)





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

# proportion true effects closer to d = 0 than d = 0.2
#  = proportion of smaller magnitude than d = 0.1
# subtract 50% to get just the positive part
prob.H0 = 2 * ( pnorm( .1 / tau ) - .5 )

# correct coarsened prior odds of H0:H1
( prior.odds = prob.H0 / (1 - prob.H0) )
# only 13%!

# where on Ioannidis' X-axis is the truth?
log( prior.odds, base = 10 )



