#% Parameters set for code
critval<-1.96
MaxIter<-10^5;
Tol<-10^(-8);
stepsize<-10^(-6);


# Height of the plot
H<-7
# Width of the plot
W<-12


# Bounds for constrained optimization for theta for correction plot
# No longer needed
if (FALSE) {
  low.bound<-(-12)
  if (application %in% c(1,2)) {
    low.bound<--100
  }
  if (application %in% c(12,3,4,5))  {
    low.bound<--2
  }
  
  u.bound<-20
}



if (identificationapproach==1) {
 # %.*sign(XX[,1]);
  Z1=Z[,1];
} else {
  Z1=sort(X/sigma);
  
}
alpha<-0.05
if (symmetric == 1) {
  Psihat_use = c(0, Psihat);
} else {
  
  if (symmetric_p == 1) {
    Psihat_use = c(Psihat[1], Psihat[2], c(1,Psihat[-c(1,2)], rev(Psihat[-c(1,2)])))
  } else {
    Psihat_use=Psihat;
  }
  
  
}
Z1_M<-rep(0,length(Z1))
Z1_L<-Z1_M
Z1_U<-Z1_M

for (n in 1:length(Z1)) {
  source("Step_function_normal_cdf.R")
  

  g_U <-function(lambda) {
    
    # MM NOTE
    # this uses the following arguments for the step function CDF:
    # X = z-score
    # cutoffs = c(1.64, 1.96)
    # theta = lambda
    # sigma = 1
    # betap = c(0.02, 0.29, 1.0)
    # symmetric = 1
    
    A<-(alpha/2 - Step_function_normal_cdf(Z1[n],lambda,1,c(Psihat_use[-c(1,2)], 1),cutoffs,symmetric))^2;
  }
  g_L <-function(lambda) {
    A<-(1-alpha/2 - Step_function_normal_cdf(Z1[n],lambda,1,c(Psihat_use[-c(1,2)], 1),cutoffs,symmetric))^2;
  }
  g_M <-function(lambda) {
    A<-(1/2 - Step_function_normal_cdf(Z1[n],lambda,1,c(Psihat_use[-c(1,2)], 1),cutoffs,symmetric))^2;
  }
  
  if (n>1) {
  Z1_U[n]<-nlm(g_U,Z1_U[n-1],gradtol=10^(-14),
               steptol=10^(-13),
               iterlim=200000000)$estimate
  Z1_L[n]<-nlm(g_L,Z1_L[n-1],gradtol=10^(-14),
               steptol=10^(-13),
               iterlim=200000000)$estimate
  Z1_M[n]<-nlm(g_M,Z1_M[n-1],gradtol=10^(-14),
                steptol=10^(-13),
                iterlim=200000000)$estimate
  }
  else {
    Z1_U[n]<-nlm(g_U,Z1[n],gradtol=10^(-14),
                 steptol=10^(-13),
                 iterlim=200000000)$estimate
    Z1_L[n]<-nlm(g_L,Z1[n],gradtol=10^(-14),
                 steptol=10^(-13),
                 iterlim=200000000)$estimate
    Z1_M[n]<-nlm(g_M,Z1[n],gradtol=10^(-14),
                 steptol=10^(-13),
                 iterlim=200000000)$estimate
  }
}


a = sapply(seq(-1,0,17),g_L)
count=1;
store<-rep(0,length(seq(-10,10,0.01)))
for (lambda in seq(-10,10,0.01)) {
  store[count]<-g_M(lambda);
  count = count+1;
}

# %% Plot original and adjusted confidence sets


Rl<-min(min(Z1_L), min(Z1)-2)-0.5;
Ru<-max(max(Z1_U), max(Z1)+2)+0.5;
R<-Ru-Rl;
H<-3.5
W<-6
library(reshape)

eps<-0.10

dat<-data.frame(original=Z1,adjusted=Z1_M,
                ytick = seq(1/n,H,H/n),
                count = 1:(length(Z1)))
df.m = melt(dat, id.vars =c("count","ytick"), measure.vars = c("original","adjusted"))
df.m$ytick[df.m$variable=="adjusted"]<-df.m$ytick[df.m$variable=="adjusted"]-eps
df.m$left<-c(Z1-critval,Z1_L)
df.m$right<-c(Z1+critval,Z1_U)

# Figure out what to plot on y axis
# default option
if (is.null(Studynames)){
  xbreaks<-1
}

xbreaks<-df.m$ytick[1:length(Z1)]

p<-ggplot(df.m, aes(ytick, value, shape = variable,colour = variable)) + 
  geom_hline(yintercept=critval,color='grey')+
  geom_hline(yintercept=0,color='grey')+
  ylim(c(min(Z1-critval,Z1_L),max(Z1+critval,Z1_U) ))+
  geom_point()+
  theme(legend.position="right", legend.direction="vertical",
        legend.title=element_blank(),
        legend.background = element_rect(linetype="solid",colour="black")) +
  scale_fill_discrete("")+
  geom_segment(aes(y=left,yend=right,x=ytick,xend=ytick,colour=variable))+
  scale_y_continuous(breaks=seq(from=floor(min(Z1-critval,Z1_L)),
                                to = ceiling(max(Z1+critval,Z1_U) ),
                                by = 2
                                ))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank())+
  scale_x_continuous(breaks=xbreaks,
                     labels =Studynames )+
  theme(legend.background = element_rect(colour = "black"),
        legend.position = c(0.85,0.5),
        axis.text.y = element_text(size=3))+
  scale_color_manual(values = c("blue", "black"))+
  coord_flip()


#add_axis(p,"x", orient = "top") 
#ggdraw(switch_axis_position(p, axis = 'x'))
#ggdraw(axis_position(p, axis = 'x',))

if (min(Z1-critval,Z1_L)<(-critval) ) {
  
  p<-p+
    geom_hline(yintercept=-critval,color='grey')
}



filepath<-paste0(pathname,'/FiguresandTables/',name, 'OriginalAndAdjusted.pdf')
save_plot(filepath,p,base_height=H,base_width=W)
filepath<-paste0(pathname,'/FiguresandTables/',name, 'OriginalAndAdjustedSlides.pdf')
save_plot(filepath,p,base_height=H,base_width=W)
# move axis to top 
# rename y axis

#%% Plot original, adjusted, and replication point estimates

if (identificationapproach==1) {
  Z2= Z[,2];
  Z2_rescaled = Z2;
  H = 3.5;
  W=6;
  
  # create dataset
  dat<-data.frame(original=Z1,adjusted=Z1_M,
                  replication = Z2_rescaled,
                  ytick = seq(1/n,H,H/n),
                  count = 1:(length(Z1)))
  # reshape into long
  df.m = melt(dat, id.vars =c("count","ytick"), measure.vars = c("original","adjusted","replication"))
  
  
  p<-ggplot(df.m, aes(ytick, value, colour = variable)) + 
    geom_line()+
    theme(legend.position="right", legend.direction="vertical",
          legend.title=element_blank(),
          legend.background = element_rect(linetype="solid",colour="black")) +
    scale_fill_discrete("")+
    scale_y_continuous(breaks=seq(from=floor(min(Z1-critval,Z1_L)),
                                  to = ceiling(max(Z1+critval,Z1_U) ),
                                  by = 2
    ))+
    theme(axis.title.x=element_blank())+
    theme(axis.title.y=element_blank())+
    scale_x_continuous(breaks=xbreaks,
                       labels = Studynames)+
    theme(legend.background = element_rect(colour = "black"),
          legend.position = c(0.85,0.5),
          axis.text.y = element_text(size=3))+
    scale_color_manual(values = c("blue", "black","grey"))+
    coord_flip()
  
  filepath<-paste0(pathname,'/FiguresandTables/',name, 'OriginalReplicationAndAdjusted.pdf')
  save_plot(filepath,p,base_height=H,base_width=W)
  filepath<-paste0(pathname,'/FiguresandTables/',name, 'OriginalReplicationAndAdjustedSlides.pdf')
  save_plot(filepath,p,base_height=H,base_width=W)
  
}

#%%plot corrected estimators and confidence sets

alpha = 0.05;


if (symmetric == 1|| symmetric_p==1) {
  xgrid<-seq(0,5,0.01)
} else {
  xgrid<-seq(-5,5,0.01)
}



# Bounds for constrained optimization for theta for correction plot
# No longer needed

if (FALSE) {
  low.bound<--2
  if (application %in% c(7,8,10,11,12,13)) {
    low.bound<--7
  }
  
  if (application %in% c(6)) {
    low.bound<--7
  }
  u.bound<-7
  
}

Theta_U_store<-rep(0,length(grid));
Theta_L_store<-Theta_U_store
Theta_M_store<-Theta_U_store
  
  
for (n in 1:length(xgrid)) {
    XX = xgrid[n];
    
    source("Step_function_normal_cdf.R")
    g_U <-function(lambda) {
      A<-(alpha/2-Step_function_normal_cdf(XX,lambda,1,c(Psihat_use[-c(1,2)], 1),cutoffs,symmetric))^2;
    }
    g_L <-function(lambda) {
      A<-(1-alpha/2-Step_function_normal_cdf(XX,lambda,1,c(Psihat_use[-c(1,2)], 1),cutoffs,symmetric))^2;
    }
    g_M <-function(lambda) {
      A<-(1/2-Step_function_normal_cdf(XX,lambda,1,c(Psihat_use[-c(1,2)], 1),cutoffs,symmetric))^2;
    }
    if (n>1) {
      Theta_U_store[n]<-nlm(g_U,Theta_U_store[n-1],gradtol=10^(-14),
                   steptol=10^(-13),
                   iterlim=200000000)$estimate
      Theta_L_store[n]<-nlm(g_L,Theta_L_store[n-1],gradtol=10^(-14),
                   steptol=10^(-13),
                   iterlim=200000000)$estimate
      Theta_M_store[n]<-nlm(g_M,Theta_M_store[n-1],gradtol=10^(-14),
                   steptol=10^(-13),
                   iterlim=200000000)$estimate
    }
    else {
      Theta_U_store[n]<-nlm(g_U,XX,gradtol=10^(-14),
                   steptol=10^(-13),
                   iterlim=200000000)$estimate
      Theta_L_store[n]<-nlm(g_L,XX,gradtol=10^(-14),
                   steptol=10^(-13),
                   iterlim=200000000)$estimate
      Theta_M_store[n]<-nlm(g_M,XX,gradtol=10^(-14),
                   steptol=10^(-13),
                   iterlim=200000000)$estimate
    }
 
}
# create data set
dat<-dat<-data.frame(U=Theta_U_store,M=Theta_M_store,
                     L = Theta_L_store,
                     xgrid=xgrid,
                     count = 1:(length(Theta_U_store)))
# reshape into long
df.m = melt(dat, id.vars =c("count","xgrid"), measure.vars = c("U","M","L"))
df.m$variable2<-factor(c(rep(1,length(Theta_U_store)),rep(2,length(Theta_U_store)),rep(1,length(Theta_U_store))
                       ))

  
  p<-ggplot(df.m,aes(x=xgrid, y = value,group=variable, size = variable2,colour ="blue"))+
      geom_abline(intercept = -critval, slope = 1, color = 'grey')+
      geom_abline(intercept = critval, slope = 1, color = 'grey')+
      geom_abline(intercept = 0, slope = 1, color = 'grey')+
      geom_line()+
      xlab('X')+
      ylab('Estimated Theta')+
      scale_x_continuous( expand = c(0,0), limits = c(min(c(0,xgrid)),ceiling(max(xgrid))),
                          breaks = seq(from = min(xgrid),
                                   to = ceiling(max(xgrid)),
                                   by =1))+
      scale_y_continuous( expand = c(0,0), limits = c(min(c(0,xgrid))-2,ceiling(max(xgrid)+2)),
                          breaks = seq(from = min(xgrid)-2,
                                     to = ceiling(max(Theta_U_store,xgrid+2)),
                                     by =1))+
      theme(legend.position=c(0.5,0.85), legend.direction="vertical",
            legend.title=element_blank(),
            legend.background = element_rect(linetype="solid",colour="grey"),
            legend.text = element_text(size = 9),
            legend.key.height = unit(0.35, "cm"))+
      scale_size_manual(values=c(0.5, 1.5),labels = c("95% Confidence Bounds", "Median Unbiased Estimator"))+
      scale_color_manual(values = c("blue", "blue"), guide=FALSE)

  
  filepath<-paste0(pathname,'/FiguresandTables/',name, 'Correction_plot.pdf')
  save_plot(filepath,p,base_height=4,base_width=4)