
############################## HELPER FUNCTIONS ##############################


sim_Z = function( .n, .mu, .tau.tilde,
                  .bp0 = 1, .bp1, .bp2 ) {
  # draw true effect vector
  # these parameters should be on Fisher's correlation scale based
  #  on pg. 60 of AK supplement
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


sim_X = function( .n, .mu, .tau, .SEs,
                  .bp0 = 1, .bp1, .bp2 ) {

  # draw true effect vector
  # these parameters should be on Fisher's correlation scale based
  #  on pg. 60 of AK supplement
  theta = rnorm( n = .n, mean = .mu,  sd = .tau )
  
  # sample from SEs
  SE = sample( .SEs, size = .n, replace = TRUE )
  
  # draw latent study for each
  Xi = rnorm( n = .n, mean = theta, sd = SE )
  
  # get Z-scores
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


# value: the value at which to evaluate 1 - ECDF
# numbers: the numbers whose inverse ECDF is being computed
inv_ecdf = function( value, numbers ) {
  # this is a function, not a number
  my.ecdf = ecdf(numbers)
  return( 1 - my.ecdf(value) )
}



# .scale: simulate on X scale vs. Z scale
# .plot.type: "scatter", "ECDF.theta", "ECDF.Z"
credibility = function( .n = 100000, .plot.n = 1000,
                        .bp0, .bp1, .bp2 = 1, .mu,
                        .tau = NA, .SEs,  # only needed if .scale = "X"
                        .tau.tilde = NA, # only needed if .scale = "Z"
                        .thresh = NA,
                        .incl.ref = TRUE,
                        .scale,
                        .plot.type ) {
  
  #browser()

  
  # TEST ONLY
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
  
  
  ##### Simulate Data - With Selection #####
  
  if ( .scale == "Z" ) {
    d = sim_Z( .n = .n, .mu = .mu, .tau.tilde = .tau.tilde,
               .bp0 = .bp0, .bp1 = .bp1, .bp2 = .bp2 )
  } else {
    d = sim_X( .n = .n, .mu = .mu, .tau = .tau, .SEs = .SEs,
               .bp0 = .bp0, .bp1 = .bp1, .bp2 = .bp2 )
  }
  
  dp = d[ d$publish == TRUE, ]
  
  ##### Simulate Data - Without Selection #####
  if ( .incl.ref == TRUE ) {
    
    if ( .scale == "Z" ) {
      D = sim_Z( .n = .n, .mu = .mu, .tau.tilde = .tau.tilde,
                 .bp0 = 1, .bp1 = 1, .bp2 = 1 )
    } else {
      D = sim_X( .n = .n, .mu = .mu, .tau = .tau, .SEs = .SEs,
                 .bp0 = 1, .bp1 = 1, .bp2 = 1 )
    }
    
    Dp = D[ D$publish == TRUE, ]  # all are published, so doesn't change dataframe at all
  }
  
  
  ##### Credibility: Proportion Above Each Threhsold #####
  # credibility for each threshold
  prop.above = vapply( X = .thresh,
                       FUN = function(t) sum( d$theta[ d$signif.pos == TRUE ] > t ) / length( d$theta[ d$signif.pos == TRUE ] ),
                       FUN.VALUE = 0.5)
  
  print( cbind( .thresh, prop.above ) )
  
  
  ##### Scatterplot #####
  # draw random sample for plotting
  if ( .plot.type == "scatter" ) {
    samp = d[ sample( 1:nrow(d), replace = TRUE, size = .plot.n ), ]
    
    library(ggplot2)
    colors = c("grey", "black")
    scatter = ggplot( data = samp, aes( x = samp$Zi, y = samp$theta, color = (samp$publish == 1) ) ) +
      geom_point( size = 1.5 ) +
      scale_color_manual(values=colors) +
      scale_x_continuous(limits=c(-6,6), breaks=seq(-6, 6, 1)) +
      scale_y_continuous(limits=c(-6,6), breaks=seq(-6, 6, 1)) +
      #geom_ribbon( aes( xmin = 1.96, ymax = Inf ), fill = "orange", alpha=0.3 ) +
      #geom_rect( aes( xmin = 1.96, xmax = Inf, ymin = -Inf, ymax = Inf ), fill = "orange", alpha=0.3 ) +
      annotate("rect", xmin=1.96, xmax=Inf, ymin=-Inf, ymax=Inf, alpha=0.2, fill="orange") +
      theme_bw() +
      xlab("Study estimate") + ylab("True effect size") + ggtitle( "Published findings only") +
      guides(color=guide_legend(title="Published"))
    
    for ( t in .thresh ) {
      scatter = scatter + geom_hline( yintercept = t, linetype=2, color = "red" )
    }
    
    plot(scatter)
  }
  
  
  ##### ECDF Of True Thetas Among Published Studies #####
  
  if ( .plot.type == "ECDF.theta" ) {
  
    colors = c("grey", "orange")
    
    E1 = ggplot( data = dp[ dp$signif.pos == TRUE, ], aes( x = theta ) ) +
      stat_function( fun = function(x) inv_ecdf( value = x, numbers = dp$theta[ dp$signif.pos == TRUE ] ),
                     aes( color = "With selectivity (MLEs)" ), lwd=1.5 ) +
      stat_function( fun = function(x) inv_ecdf( value = x, numbers = Dp$theta[ Dp$signif.pos == TRUE ] ),
                     aes( color = "No selectivity (publish everything)" ), lwd=1 ) +
      theme_bw() +
      xlab("True effect size") + ylab("Proportion above") + ggtitle( "1 - ECDF") +
      scale_x_continuous(limits=c(-6,6), breaks=seq(-6, 6, 1)) +
      scale_y_continuous(limits=c(0,1), breaks=seq(0, 1, .05)) +
      scale_color_manual( name=" ", values=colors )
    
    for ( t in .thresh ) {
      E1 = E1 + geom_vline( xintercept = t, linetype=2, color = "red" )
    }
    plot(E1)
  }
  
  
  
  ##### ECDF Of Z-Scores #####
  
  if ( .plot.type == "ECDF.Z" ) {
    E2 = ggplot( data = dp, aes( dp$Zi ) ) +
      stat_ecdf(geom = "step") +
      stat_ecdf( data = d, aes(x = d$Zi), color = "blue" ) +  # ECDF without selectivity (all findings)
      theme_bw() +
      scale_color_manual(values=colors) +
      xlab("Z-score") + ylab("Proportion below")
    plot(E2)
  }
  
  # return dataset invisibly
  invisible(d)
  
  
  ##### Sanity checks #####
  
  # how often are they significant?
  #prop.table( table( Zi.star.abs > 1.96 ) )  # 8% of true ones are significant
  
  # prop.table( table( abs(Zi) > 1.96 ) )  # 63% of published ones are significant
  # prop.table( table( s$t > 1.96 ) )  # cf. 63% in data
  
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
  
  
  # now plot just the published, significant, positive ones
  # allows us to eyeball our credibility conditional prob:
  #  number of black points / total points on this graph
  # plot( d2$Zi.star, d2$theta, xlim = c(1.8,8), ylim = c(-4,8), col = ifelse(d2$theta < 0, "red", "black") )
  # abline(a=0, b=0, col="red", lty=2)
  
}

