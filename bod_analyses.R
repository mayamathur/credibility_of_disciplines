
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



############################## SANITY CHECK: SIMULATION ##############################

######## Parameters (Used Throughout) ########

# from running AK with Z-scores
bp0 = 1  # prob of publishing significant result (baseline)
bp1 = 0.3864854  
bp2 = 0.4464266
mu = 0
tau.tilde = 2.118972  

# from running AK with raw effect sizes
bp0 = 1
bp1 = 0.441
bp2 = 0.491
tau = 0.584
mu = 0


# BOOKMARK: WAS GOING TO TRY RUNNING CREDIBILITY THING USING NEW HELPER CODE
# THEN READ IOANNIDIS' PAPER TO UNDERSTAND WHAT EXACTLY HE TRIED TO ESTIMATE AND WHY IT'S WRONG

source("helper.R")
credibility( .n = 5000,
               .plot.n = 5000,
               .bp0 = bp0,
               .bp1 = bp1,
               .bp2 = 1,
               .mu = mu,
               .SEs = s$SE,
               .tau = tau,
               .thresh = c(0, .1, .5),
               .incl.ref = TRUE,
               .scale = "X",
               .plot.type = "ECDF.Z" )








######## Sanity Check: Simulate Using MLEs ########


######## The Answer ########
# P( theta > 0 | Z.orig > 1.96 )
prop.table( table( d2$theta > 0) )

# shouldn't this be equivalent?
#sum( theta[ Zi > 1.96 ] > 0 ) / length( theta[ Zi > 1.96 ] )
# THIS ONE IS WRONG BECAUSE ZI IS DIFFERNT LENGTH, SO GETS RECYCLED

# equivalent but shorter:
sum( theta[ Zi.star > 1.96 ] > 1 ) / length( theta[ Zi.star > 1.96 ] )

# look at different threshold to see if it depends on selectivity
# yes - from 91% (no selectivity) to 99% (horrible selectivity)
sum( theta[ Zi.star > .2 & publish == TRUE ] > 0 ) / length( theta[ Zi.star > .2 & publish == TRUE ] )

# BOOKMARK
# SHOW TYLER





######## Sanity Check: Compare Entire Theoretical CDF to ECDF ########

##### Fn: Given parameters, return theoretical CDF values ######

get_TCDF = function( Z.grid, .betap, .cutoffs, .tau ) {
  
  integrand = function( theta.star, Z ) {
    # P( Z-score < Z | Theta = theta* )
    term1 = Step_function_normal_cdf( X = Z, theta = theta.star, sigma = 1,
                                      betap = .betap, cutoffs = .cutoffs,
                                      symmetric = 1 )
    
    # f( theta* )
    term2 = (1 / .tau) * dnorm( theta.star / .tau )
    
    # P( Di = 1 | Z )
    #if ( abs(Z) )
    
    return( term1 * term2 )
  }

  ##### Compute TCDF values ######
  TCDF = c()
  for (i in 1:length(Z.grid)) {
    TCDF[i] = integrate( f = function(x) integrand( theta.star = x, Z = Z[i] ),
                         lower = -Inf, upper = Inf )$value
  }
  
  # TCDF = vapply( X = Z.grid, FUN = integrate( f = function(x) integrand( theta.star = x, Z = Z[i] ),
  #                                      lower = -Inf, upper = Inf )$value )
  return(TCDF)
}


# ##### MLE ######
# Z.grid = seq(-10,10,.1)
# TCDF.MLE3 = get_TCDF( Z.grid = Z.grid, .betap = c( bp1, bp2, 1), .cutoffs = c(1.64, 1.96), .tau = tau.tilde )
# 
# 
# ##### MLE, 5 cutoffs ######
# Z.grid = seq(-10,10,.1)
# TCDF.MLE5 = get_TCDF( Z.grid = Z.grid, .betap = c( 0.0705629, 0.09186333, 0.225171, 0.3152829, 0.5089435, 1), .cutoffs = c(1.64, 1.96, 3, 4, 5), .tau = 1.598206 )



##### No-selection #####
TCDF.NS = get_TCDF( Z.grid = Z.grid, .betap = c( 1, 1, 1), .cutoffs = c(1.64, 1.96), .tau = tau.tilde )


# ##### Horrible-selection version #####
# TCDF.HS = get_TCDF( Z.grid = Z.grid, .betap = c( 0, 0, 1), .cutoffs = c(1.64, 1.96), .tau = tau.tilde )



# plot 4 ECDFs together
plot( ecdf( s$tr ), xlim=c(-10,10), lwd = 5 )
abline(v=c(-1.96, -1.64, 1.64, 1.96), lty=2, col="red" )
# lines( Z, TCDF.MLE, col = "red", lwd=2 )  # MLE integral
# lines( Z, TCDF.MLE5, col = "red", lwd=2 )  # MLE integral with 5 cutoffs doesn't help
lines( ecdf( Zi ), xlim=c(-10,10), col = "green", lwd=2 )  # MLE simulation
lines( Z, TCDF.NS, xlim=c(-10,10), col = "blue", lwd=2 )  # no-selection integral
#lines( Z, TCDF.HS, xlim=c(-10,10), col = "blue", lwd=2 )  # horrible-selection integral




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                               DEBUGGING 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# find out why simulation and CDF don't match 

# bp1 = 1
# bp2 = 1

bp0 = 1  # baseline scale factor - doesn't change anything
bp1 = 0.3864854  
bp2 = 0.4464266

# bp1 = .1
# bp2 = .2


######## 1. Simulate at specific value of theta for P( Z < 1.96 | Theta = theta )

# number of studies to draw (arbitrary)
n = 500000

# draw true effect vector
# these parameters should indeed be on Fisher's correlation scale based
#  on pg. 60 of AK supplement

theta = rnorm( n = n, mean = mu,  sd = tau.tilde )
#theta = 1.96

# draw latent study for each
Zi.star = rnorm( n = n, mean = theta, sd = 1 )

# draw publication indicator
library(car)

publish.prob = rep( bp0, n )
publish.prob[ abs(Zi.star) < 1.64 ] = bp1 * bp0
publish.prob[ ( abs(Zi.star) >= 1.64 ) & ( abs(Zi.star) < 1.96 ) ] = bp2 * bp0
publish = rbinom( size = 1, n = n, prob = publish.prob )

aggregate( abs(Zi.star) ~ publish.prob, FUN = mean)
aggregate( abs(Zi.star) ~ publish.prob, FUN = max)
aggregate( abs(Zi.star) ~ publish.prob, FUN = min)

Zi = Zi.star[publish==TRUE]


# P( Z < 1.96 | Theta = theta )
sum( Zi <= 1.96 ) / length(Zi)
# for theta = 0, should be 1 - 0.025 = 0.975: matches




######## 1. Compare to theory
Step_function_normal_cdf( X = 1.96, theta = theta, sigma = 1,
                          betap = c(bp1, bp2, 1), cutoffs = c(1.64, 1.96),
                          symmetric = 1 )

# ~~~~~~~ AGREES FOR VARIOUS CHOICES OF THETA AND BPS
# even theta = 1.96, a potentially problematic area




######## 2. Simulation (integrate over thetas)

######## 
# now try integrating over thetas
integrand = function( theta.star ) {
  # P( Z-score < Z | Theta = theta* )
  term1 = Step_function_normal_cdf( X = 1.96, theta = theta.star, sigma = 1,
                                    betap = c(bp1, bp2, 1), cutoffs = c(1.64, 1.96),
                                    symmetric = 1 )
  
  # f( theta* )
  term2 = (1 / tau.tilde) * dnorm( theta.star / tau.tilde )
  
  # P( Di = 1 | Z < 1.96 )
  
  return( term1 * term2 )
}

integrate( f = function(x) integrand( theta.star = x ),
                       lower = -Inf, upper = Inf )$value

# ~~~~~~~~~~~~ DISAGREES!!!!!!!!!!!



######## 3. Do I have distribution of theta right?

threshold = .39
sum( theta < threshold ) / length(theta)

integrate( f = function(x) (1 / tau.tilde) * dnorm( x / tau.tilde ),
           lower = -Inf, upper = threshold )$value
# AGREES



######## 4. Do my own "numerical integration" using simulated thetas

# thetas for published Z's
theta.pub = theta[ publish == TRUE ]

cond.prob = c()

for ( i in 1:length(theta.pub) ) {
  
  # compute conditional prob at this theta
  cond.prob[i] = Step_function_normal_cdf( X = 1.96, theta = theta.pub[i], sigma = 1,
                            betap = c(bp1, bp2, 1), cutoffs = c(1.64, 1.96),
                            symmetric = 1 )
}

# average them
mean( cond.prob )

# THIS AGREES!
# PROBLEM: INTEGRAL IS OVER MARGINAL DIST OF THETAS, BUT I THINK IT HAS TO BE OVER 
#  DIST OF THETAS FOR PUBLISHED STUDIES ONLY. 
#  THIS IS WHY THE INTEGRAL AND SIM AGREE EXACTLY WHEN THERE'S NO SELECTIVITY: 
#  THEN THE TWO THETA DISTS ARE THE SAME. 





# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                               CREDIBILITY: P( THETA > 0 | Z > 1.96 )
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


######## Numerator Conditional Probability ########

# function to integrate to get P( Z.orig > 1.96 | Theta > 0)

integrand = function( theta.star ) {
  
  # P( Z.orig > 1.96 | Theta = theta* )
  term1 = 1 - Step_function_normal_cdf( X = 1.96, theta = theta.star, sigma = 1,
                                        #betap = c(1, 1, 1), cutoffs = c( 1.64, 1.96 ),
                                        betap = c(bp1, bp2, 1), cutoffs = c( 1.64, 1.96 ),
                                        symmetric = 1 )
  
  
  # f( theta* | Theta > 0 )
  # note presence of indicator for whether theta.star > 0
  # hence why we can integrate over whole space of Theta
  tau = tau.tilde  # NOW USING THE TAU FOR THE Z-SCORE; SEE MY NOTE IN STEP_FUNCTION ABOUT THIS
  term2 = 2 * (1 / tau) * dnorm( theta.star / tau ) * (theta.star > 0)
  # this is a half-normal distribution (i.e., folded normal with mean=0)
  
  return( term1 * term2 )
}

( p1 = integrate( f = function(x) integrand( theta.star = x ), lower = -Inf, upper = Inf )$value )


# as sanity check, changing to no selectivity in integrand 
#  gives p1 = 0.04
#  and p.temp = 0.01 below
#  so we exactly recover the alpha level by summing the two ways to reject, as expected




######## Numerator Marginal Probability ########

# P( theta > 0 )
p2 = 0.5  # because N(0, tau)



######## Denominator Conditional Probability ########

# P( Z.orig > 1.96 | Theta < 0 )
integrand3 = function( theta.star ) {
  
  # P( Z-score > 0 | Theta = theta* )
  term1 = 1 - Step_function_normal_cdf( X = 1.96, theta = theta.star, sigma = 1,
                                        #betap = c(1, 1, 1), cutoffs = c( 1.64, 1.96 ),
                                        betap = c(bp1, bp2, 1), cutoffs = c( 1.64, 1.96 ),
                                        symmetric = 1 )
  
  
  # f( theta* )
  tau = tau.tilde
  term2 = (1 / tau) * dnorm( theta.star / tau )
  
  return( term1 * term2 )
}

( p.denom = integrate( f = function(x) integrand3( theta.star = x ),
                       lower = -Inf, upper = Inf )$value )



####### Sanity check: does 2 * p.denom match empirical? ########
2 * p.denom
sum( s$t > 1.96 ) / nrow(s)


###### The Answer #######
# P( theta > 0 | Z.orig > 1.96 )
( p1 * p2 ) / p.denom
# 98% 



##### compare to simulation #####
sum( theta[ Zi > 1.96 ] > 0 ) / length( theta[ Zi > 1.96 ] )



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                               P( THETA > THRESHOLD | Z > 1.96 )
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# BOOKMARK
# DOESN'T SEEM TO CHANGE BASED ON BP PARAMETERS


# from passing raw effect sizes and SEs to AK code
tau = 0.528
bp1 = 0.357
bp2 = 0.426
threshold = 0.2  # on Cohen's d scale

# # no selection
# bp1 = 1
# bp2 = 1

# # horrible selection
# bp1 = 0
# bp2 = 0

# number of studies to draw (arbitrary)
n = 100000


# draw true effect vector
# these parameters should indeed be on Fisher's correlation scale based
#  on pg. 60 of AK supplement
theta = rnorm( n = n, mean = mu,  sd = tau )

# draw SE for latent study
# AK included the interactions, so let's start there
# get vector of SEs from RPP
SEs = sample( s$SE, size = n, replace = TRUE )

# draw latent study for each
Xi.star = rnorm( n = n, mean = theta, sd = SEs )

# latent study estimates vs. their true thetas
# no selection yet
# plot( theta, Xi.star, xlim = c(-2,2), ylim = c(-2,2) )
# abline(a=0, b=1, col="red", lty=2)



# get Z-stats
Zi.star.abs = abs( Xi.star / SEs )
# how often are true one significant?
prop.table( table( Zi.star.abs > 1.96 ) )

# draw publication indicator
library(car)
publish.prob = rep( bp0, n )
publish.prob[ Zi.star.abs < 1.64 ] = bp1 * bp0
publish.prob[ ( Zi.star.abs >= 1.64 ) & ( Zi.star.abs < 1.96 ) ] = bp2 * bp0
publish = rbinom( size = 1, n = n, prob = publish.prob )

# sanity check
Zi = Xi.star[publish==TRUE] / SEs[publish==TRUE]
prop.table( table( abs(Zi) > 1.96 ) )  # 64% of published ones are significant
# 2 * prop.table( table( Zi > 1.96 ) )  # should be similar to above by symmetry
# empirical (64%) is close to the above 63%
prop.table( table( s$t > 1.96 ) )

# sanity check
aggregate( publish ~ publish.prob, FUN = mean )
aggregate( abs(Xi.star) ~ publish.prob, FUN = mean )
aggregate( abs(Xi.star) ~ publish, FUN = mean )  # pub bias in effect sizes

# ones that survive to publication
Xi = Xi.star[ publish == 1 ]
length(Xi) / n  # proportion that survived = 62%
mean(abs(Xi)); mean(abs(Xi.star)) # as expected, the average estimated effect size is larger in ones that survive

mean( abs( theta[publish == TRUE] ) ); mean( abs( theta ) ) # as are the true effect sizes among published ones


# sanity check: make plot similar to Fig 4 (pg 23)
# under selection, this will have discontinuities
# plot( Xi, theta[publish == TRUE], xlim = c(-1,1), ylim = c(-1,1), col = ifelse( Zi > 1.96, "red", "black") )
# abline(a=0, b=1, col="red", lty=2)


# put things in dataframe
d = data.frame( Xi.star, publish, SEs, theta )
d$keep = ( d$Xi.star / d$SEs ) > 1.96  # is study positive and significant?

# change this to assess effect sizes other than 0
d$true.pos = d$theta > threshold

# keep only published, significant, positive-estimate studies
# last condition technically redundant since all significant ones are published
d2 = d[ d$keep == 1 & d$publish == 1, ]
nrow(d2) / n  # proportion kept: 16% regardless of selectivity (makes sense)

# now plot just the published, significant, positive ones
# allows us to eyeball our credibility conditional prob:
#  number of black points / total points on this graph
# plot( d2$Xi.star, d2$theta, xlim = c(0,1), ylim = c(-1,1), col = ifelse(d2$theta < threshold, "red", "black") )
# abline(a=threshold, b=0, col="red", lty=2)


######## The Answer ########
# P( theta > 0 | Z.orig > 1.96 )
prop.table( table( d2$theta > threshold) )
# 64% for d = 0.5
# 97% for d = 0.1

# bookmark: doesn't depend on selection params....
# HYPOTHESIS: THIS IS BECAUSE NO SELECTION OPERATES ON THE Z > 1.96 STUDIES
# SINCE THEY ARE ALL PUBLISHED. TRY USING DIFFERENT THRESHOLD TO CHECK THIS. 


