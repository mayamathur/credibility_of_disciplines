
############################## HELPER FUNCTIONS ##############################


######## Fns: Ioannidis FRPs ########
# see Ioannidis' first equation
frp_0 = function( O, alpha, pwr ) {
  ( O * alpha ) / ( O * alpha + pwr )
}

# see Ioannidis pg 7
frp_thresh = function( O, alpha, pwr.s, prob.s, pwr.l, prob.l ) {
  ( O * alpha + pwr.s * prob.s ) / ( O * alpha + pwr.s * prob.s + pwr.l * prob.l )
}



######## Fn: Simulate Z-Score Effect Sizes ########

# .bp0 = baseline publication probability (for |Z| > 1.96)
# .bp1 = relative publication probability for nonsignificant Z
# .bp2 = relative publication probability for marginally significant Z

sim_Z = function( .n, .mu, .tau.tilde,
                  .bp0 = 1, .bp1, .bp2 ) {
  
  # draw true effect vector on Z-score scale
  theta = rnorm( n = .n, mean = .mu,  sd = .tau.tilde )
  
  # draw latent study for each
  Zi = rnorm( n = .n, mean = theta, sd = 1 )
  
  d = data.frame( theta, Zi )
  
  # draw publication indicator
  library(car)
  d$publish.prob = rep( .bp0, .n )
  d$publish.prob[ abs( d$Zi ) < 1.64 ] = .bp1 * .bp0
  d$publish.prob[ ( abs( d$Zi ) >= 1.64 ) & ( abs( d$Zi ) < 1.96 ) ] = .bp2 * .bp0
  d$publish = rbinom( size = 1, n = .n, prob = d$publish.prob )
  
  d$signif.pos = d$Zi > 1.96
  
  invisible( return(d) )
}


######## Fn: Simulate Raw Effect Sizes ########

# if using AK's parameter estimates on RPP data, 
#  parameters are on Fisher's correlation scale based
#  on pg. 60 of AK supplement

sim_X = function( .n, .mu, .tau, .SEs,
                  .bp0 = 1, .bp1, .bp2 ) {

  # draw true effect vector
  theta = rnorm( n = .n, mean = .mu,  sd = .tau )
  
  # sample from SEs
  SE = sample( .SEs, size = .n, replace = TRUE )
  
  # draw latent study for each
  Xi = rnorm( n = .n, mean = theta, sd = SE )
  Zi = Xi / SE
  
  d = data.frame( theta, Xi, SE, Zi )
  
  # draw publication indicator
  library(car)
  d$publish.prob = rep( .bp0, .n )
  d$publish.prob[ abs( d$Zi ) < 1.64 ] = .bp1 * .bp0
  d$publish.prob[ ( abs( d$Zi ) >= 1.64 ) & ( abs( d$Zi ) < 1.96 ) ] = .bp2 * .bp0
  d$publish = rbinom( size = 1, n = .n, prob = d$publish.prob )
  
  d$signif.pos = d$Zi > 1.96
  
  invisible( return(d) )
}


######## Fn: Evaluate 1 - ECDF at a Value ########

# value: the value at which to evaluate 1 - ECDF
# numbers: the numbers whose inverse ECDF is being computed

inv_ecdf = function( value, numbers ) {
  my.ecdf = ecdf(numbers) # this is a fn, not a scalar
  return( 1 - my.ecdf(value) )
}


# .n = sample size for simulation
# .plot.n = random subset size for plots
# .prox = proximity of reported ES to truth (as absolute deviation from truth)
# .scale = simulate on "X" scale vs. "Z" scale
# .plot.type = "scatter", "ECDF.theta", "ECDF.Z"
# .incl.ref = should reference line (no selectivity) be included in CDF plots?

credibility = function( .n = 100000, .plot.n = 1000,
                        .bp0, .bp1, .bp2 = 1, .mu,
                        .tau = NA, .SEs,  # only needed if .scale = "X"
                        .tau.tilde = NA, # only needed if .scale = "Z"
                        .thresh = NA,
                        .prox = NA, 
                        .incl.ref = TRUE,
                        .scale,
                        .plot.type ) {

  # TEST ONLY
  #  browser()
  # .n = 1000
  # .plot.n = 800
  # .bp0 = bp0
  # .bp1 = bp1
  # .bp2 = 1
  # .mu = mu
  # .tau.tilde = tau.tilde
  # .SEs = s$SE * 4
  # .tau = tau
  # .thresh = c(0, .1)
  # .incl.ref = TRUE
  # .scale = "X"
  

  ##### Simulate Data With Selection #####
  
  if ( .scale == "Z" ) {
    d = sim_Z( .n = .n, .mu = .mu, .tau.tilde = .tau.tilde,
               .bp0 = .bp0, .bp1 = .bp1, .bp2 = .bp2 )
  } else {
    d = sim_X( .n = .n, .mu = .mu, .tau = .tau, .SEs = .SEs,
               .bp0 = .bp0, .bp1 = .bp1, .bp2 = .bp2 )
  }
  
  dp = d[ d$publish == TRUE, ]
  

  ##### Proximity: Proportion Within Caliper of Truth #####
  if ( .scale == "X" & !is.na(.prox) ) {
   
    # uses caliper as proportion
    # prop.calip = vapply( X = .prox,
    #                      FUN = function(t) sum( abs( temp$Xi - temp$theta ) <= t * abs( temp$theta ) ) /
    #                        nrow(temp),
    #                      FUN.VALUE = 0.5)
    
    # uses caliper as absolute deviation
    prop.calip = vapply( X = .prox,
                         FUN = function(t) sum( abs( dp$Xi - dp$theta ) <= t ) /
                           nrow(dp),
                         FUN.VALUE = 0.5 )
    print( cbind( .prox, prop.calip ) )
  }
  
  if ( .scale == "Z" & !is.na(.prox) ) {

    prop.calip = vapply( X = .prox,
                         FUN = function(t) sum( abs( dp$Zi - dp$theta ) <= t ) /
                           nrow(dp),
                         FUN.VALUE = 0.5 )
    print( cbind( .prox, prop.calip ) )
  }

  
  ##### Avg Distance: Distance of Reported Value from Truth #####
  if ( .scale == "X" ) {
  
    # uses caliper as absolute deviation
    dist = abs( dp$Xi - dp$theta )
      
    cat( paste( "\n\nMedian distance of published finding from truth: ", round( median(dist), 2 ),
                "\n", sep = "" ) )
    cat( paste( "\n\nMedian distance (no selectivity): ", round( median( abs( d$Xi - d$theta ) ), 2 ),
                "\n\n", sep = "" ) )
  }
  
  if ( .scale == "Z" ) {

    # uses caliper as absolute deviation
    dist = abs( dp$Zi - dp$theta )
    
    cat( paste( "\n\nMedian distance of published finding from truth: ", median(dist),
                "\n\n", sep = "" ) )
  }
  
  ##### Credibility: Proportion Above Each Threhsold #####
  # credibility for each threshold
  temp = d[ d$signif.pos == TRUE, ]
  prop.above = vapply( X = .thresh,
                       FUN = function(t) sum( temp$theta > t ) / nrow(temp),
                       FUN.VALUE = 0.5 )
  
  print( cbind( .thresh, prop.above ) )
  
  
  ##### Scatterplot #####
  # draw random sample for plotting
  if ( .plot.type == "scatter" ) {
    samp = d[ sample( 1:nrow(d), replace = TRUE, size = .plot.n ), ]
    
    # add new column for plotting to ensure it's on same scale as theta
    if ( .scale == "Z" ) samp$xvar = samp$Zi
    else if ( .scale == "X" ) samp$xvar = samp$Xi
    
    library(ggplot2)
    colors = c("grey", "black")
    scatter = ggplot( data = samp, aes( x = samp$xvar,
                                        y = samp$theta, color = (samp$publish == 1) ) ) +
      geom_abline( intercept=0, slope=1, color = "grey" ) +
      geom_point( size = 1.5 ) +
      scale_color_manual(values=colors) +
      scale_x_continuous(limits=c(-6,6), breaks=seq(-6, 6, 1)) +
      scale_y_continuous(limits=c(-6,6), breaks=seq(-6, 6, 1)) +
      annotate("rect", xmin=1.96, xmax=Inf, ymin=-Inf, ymax=Inf, alpha=0.2, fill="orange") +
      theme_bw() +
      xlab( paste( "Study estimate on ", .scale, " scale", sep = "" ) ) +
      ylab( paste( "True effect size on ", .scale, " scale", sep = "" ) ) +
      ggtitle( "All studies (shaded: Z > 1.96)") +
      guides(color=guide_legend(title="Published"))
    
    for ( t in .thresh ) {
      scatter = scatter + geom_hline( yintercept = t, linetype=2, color = "red" )
    }
    
    plot(scatter)
  }
  
  
  ##### Distance-From-Truth Scatterplot #####
  # draw random sample for plotting
  if ( .plot.type == "dist" ) {
  
    samp = d[ sample( 1:nrow(d), replace = TRUE, size = .plot.n ), ]
    
   # browser()
    
    # add new column for plotting to ensure it's on same scale as theta
    if ( .scale == "Z" ) samp$xvar = samp$Zi
    else if ( .scale == "X" ) samp$xvar = samp$Xi
    
    library(ggplot2)
    colors = c("grey", "black")
    scatter = ggplot( data = samp, aes( x = samp$xvar, y = abs( samp$xvar - samp$theta ),
                                        color = (samp$publish == 1) ) ) +
      geom_point( size = 1.5 ) +
      scale_color_manual(values=colors) +
      scale_x_continuous(limits=c(-6,6), breaks=seq(-6, 6, 1)) +
      scale_y_continuous(limits=c(-6,6), breaks=seq(-6, 6, 1)) +
      annotate("rect", xmin=1.96, xmax=Inf, ymin=-Inf, ymax=Inf, alpha=0.2, fill="orange") +
      theme_bw() +
      xlab( paste( "Study estimate on ", .scale, " scale", sep = "" ) ) +
      ylab( paste( "Distance from truth on ", .scale, " scale", sep = "" ) ) +
      ggtitle( "All studies (shaded: Z > 1.96)") +
      guides(color=guide_legend(title="Published"))
    
    for ( t in .thresh ) {
      scatter = scatter + geom_hline( yintercept = t, linetype=2, color = "red" )
    }
    
    plot(scatter)
  }
  
  
  ##### ECDF Of True Thetas Among Z > 1.96 Studies #####
  
  if ( .plot.type == "ECDF.theta" ) {
  
    colors = c("grey", "orange")
    
    E1 = ggplot( data = dp[ dp$signif.pos == TRUE, ], aes( x = theta ) ) +
      
      stat_function( fun = function(x) inv_ecdf( value = x, numbers = dp$theta[ dp$signif.pos == TRUE ] ),
                     aes( color = "With selectivity (MLEs)" ), lwd=1.5 ) +
      
      stat_function( fun = function(x) inv_ecdf( value = x, numbers = Dp$theta[ Dp$signif.pos == TRUE ] ),
                     aes( color = "No selectivity (publish everything)" ), lwd=1 ) +
      theme_bw() +
      xlab("True effect size") + ylab("Proportion above") +
      ggtitle( "1 - ECDF of true effect sizes among published, significant studies") +
      scale_x_continuous(limits=c(-6,6), breaks=seq(-6, 6, 1)) +
      scale_y_continuous(limits=c(0,1), breaks=seq(0, 1, .05)) +
      scale_color_manual( name=" ", values=colors )
    
    for ( t in .thresh ) {
      E1 = E1 + geom_vline( xintercept = t, linetype=2, color = "red" )
    }
    plot(E1)
  }
  
  
  
  ##### ECDF Of All Z-Scores #####
  
  if ( .plot.type == "ECDF.Z" ) {
    
    # get Ioannidis data for empirical comparison
    setwd("~/Dropbox/Personal computer/Independent studies/BOD (believability of disciplines)/bod_git")
    s = read.csv("Data from Ioannidis/full_data.csv")
    
    colors = c("red", "grey", "orange")

    E2 = ggplot( data = dp, aes( dp$Zi ) ) +
<<<<<<< HEAD
      stat_ecdf( data = d, aes( x = d$Zi, color = "No selectivity (publish everything)" ),
                 show.legend = TRUE, lwd = 1 ) +  # ECDF without selectivity (all findings)
      
      stat_ecdf( aes( color = "With user-specified selectivity" ),
                 show.legend = TRUE, lwd = 1.5 ) +  # ECDF of simulated data
      
      stat_ecdf( data = s, aes( x = s$tr, color = "Empirical (random signs)" ),
                 show.legend = TRUE, lwd = 1 ) +  # ECDF without selectivity (all findings)
=======
      stat_ecdf( data = d, aes( x = d$Zi, color = "No selectivity (publish everything)" ), show.legend = TRUE, lwd = 1 ) +  # ECDF without selectivity (all findings)
      stat_ecdf( aes( color = "With user-specified selectivity" ), show.legend = TRUE, lwd = 1.5 ) +  # ECDF of simulated data
      stat_ecdf( data = s, aes( x = s$tr, color = "Empirical (random signs)" ), show.legend = TRUE, lwd = 1 ) +  # ECDF without selectivity (all findings)
>>>>>>> 3327ed475b08aa40e5ed050cfabb1358f6b83fdb
      theme_bw() +
      scale_color_manual( name = " ", values = colors) +
      scale_x_continuous(limits=c(-6,6), breaks=seq(-6, 6, 1)) +
      scale_y_continuous(limits=c(0,1), breaks=seq(0, 1, .05)) +
      geom_vline( xintercept = c(-1.64, 1.64, 1.96, -1.96), linetype=2, color = "red" ) +
      xlab("Z-score") + ylab("Proportion below") +
      ggtitle("Model fit: ECDF of published Z-scores")
    plot(E2)
    
    # base R plotting: runs a bit faster
    # plot( ecdf( s$tr ), xlim=c(-6,6), lwd = 5 )
    # abline(v=c(-1.96, -1.64, 1.64, 1.96), lty=2, col="red" )
    # lines( ecdf( d$Zi ), xlim=c(-6,6), col = "green", lwd=2 )  # no selection
    # lines( ecdf( dp$Zi ), xlim=c(-6,6), col = "orange", lwd=2 )  # with selection
  }
  
  # return dataset invisibly
  invisible(d)
  
  
  ##### Sanity checks #####
  
  # sanity check
  # aggregate( publish ~ publish.prob, FUN = mean )
  # aggregate( abs(Zi.star) ~ publish.prob, FUN = mean )
  # aggregate( abs(Zi.star) ~ publish, FUN = mean )  # pub bias in effect sizes
  
  # mean effect size for significant & positive vs. not
  #aggregate( theta ~ ( Zi.star > 1.96  ), FUN=mean )
  
  # ones that survive to publication
  # length(Zi) / n  # proportion that survived = 41%
  # mean(abs(Zi)); mean(abs(Zi.star)) # as expected, the average Z-score is larger in ones that survive
  # mean( abs( theta[publish == TRUE] ) ); mean( abs( theta ) ) # as are the true effect sizes among published ones
  
}




##### Convert Between r (Correlation) and z (Fisher's) Units ##### 
r_to_z = Vectorize( function(r) {
  .5 * ( log(1 + r) - log(1 - r) )
}, vectorize.args = "r" )

z_to_r = Vectorize( function(z) {
  ( exp( 2 * z ) - 1 ) / ( exp( 2 * z ) + 1 )
}, vectorize.args = "z" )

# convert Cohen's d to r
# assumes equal sample sizes in each group
# see Borenstein text
d_to_r = Vectorize( function(d) {
  d / sqrt(d^2 + 4)
}, vectorize.args = "d" )


r_to_d = Vectorize( function(r) {
  (2 * r) / sqrt(1 - r^2)
}, vectorize.args = "r" )



