

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                               LOAD IOANNIDIS DATA
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

setwd("~/Dropbox/Personal computer/Independent studies/BOD (believability of disciplines)/bod_git")

library(R.matlab)
d1 = readMat("Data from Ioannidis/data.mat")
data = d1$D
data = data[, , 1]

# look at other fields that were wrong length
data$JournalTypes
data$JournalNames
data$Comment

d = data.frame( tvalues = data$tvalues, df = data$df )

# significant?
d$t = d$tvalues
d$pval = 2 * ( 1 - pt( q = d$t, df = d$df ) )
d$rej = d$pval < 0.05


##### Probabilistic Mixture Assignment of Effect Sizes #####
# per Ioannidis paper

# Cohen's d assuming 2-sample t-test
d$D2 = ( 2 * d$tvalues ) / sqrt( d$df + 2 )
# Cohen's d assuming 1-sample t-test
d$D1 = d$tvalues / sqrt( d$df + 1 )

# probability it's a 1-sample t-test given df
d$p.one.samp = ifelse( d$df <= 10, 0.93, 0.72 )

# mixture thing to estimate effect size
d$D = d$D1 * d$p.one.samp + d$D2 * (1 - d$p.one.samp) 

# total sample size
N1 = d$df + 1
N2 = d$df + 2
d$N = d$N1 * d$p.one.samp + d$N2 * (1 - d$p.one.samp) 

# get SE for Cohen's d
# because D / SE = t
d$SE = d$D / d$t


##### Randomly Assign Signs #####
set.seed(451); sign = ifelse( rbinom( n = nrow(d), size=1, prob = 0.5 ), 1, -1 )
d$Dr = sign * d$D
d$tr = sign * d$t



############################## EXPLORE DATA ############################## 

# mean df = 67
mean(d$df)

# reproduce their median effect sizes
# reported: 0.932 and 0.237
# matches :) 
aggregate( D ~ rej, FUN = median, data = d )

# 64% were significant
prop.table( table( d$rej) )

# number with df > 30 (good for Z approximation)
table(sign(d$tvalues)) # already sign-normalized because none are negative

# histogram minus some extreme outliers
hist(d$tvalues[ d$tvalues<20 ] )
summary(d$tvalues)
prop.table(table(d$tvalues>20))  # 99% are t < 20


##### Exclusions #####
# keep non-extreme ones and ones with 
d$exclude = ( d$t > 6 ) | is.na( d$SE )
prop.table( table(d$exclude) )  # proportion exclusions
d = d[ d$exclude == FALSE, ]



############################## RANDOM SAMPLE FOR CODE DEVELOPMENT ############################## 

n = 5000
set.seed(451); s = d[ sample( 1 : nrow(d), size = n ), ]

write.csv( s, "random_sample.csv" )
write.csv( d, "full_data.csv" )
