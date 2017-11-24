
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



