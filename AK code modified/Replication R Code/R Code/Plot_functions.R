replication_plot<-function(X, xaxisname,yaxisname) {
  
  ll<-floor(min(min(X)));
  uu<-ceiling(max(max(X)));
  significant<-(abs(X[,1])>critval);
  dat<-data.frame(xvar = X[,1],
                  yvar = X[,2],
                  significant = as.factor(significant))
  dat<-dat[order(dat$significant),]
  

  
  # Plot scatter and scatterslides (with letters A,B)    
  p<- ggplot(dat, aes(x=xvar, y=yvar)) +
    geom_vline(xintercept =critval,color="grey")+
    geom_vline(xintercept =-critval,color="grey")+
    geom_hline(yintercept =critval,color="grey")+
    geom_hline(yintercept =-critval,color="grey")+
    geom_abline(intercept = 0,slope=1,color="grey")+
    xlab(xaxisname)+
    ylab(yaxisname )+
    xlim(c(ll,uu))+
    ylim(c(ll,uu))+
    geom_point(shape=shape,size = 2,aes(colour = significant,
                                        fill = significant))+ 
    scale_fill_manual(values=c("grey", "blue")) + 
    scale_colour_manual(values=c("grey", "blue"))
#scale_x_continuous(breaks = round(seq(min(dat$xvar), max(dat$xvar), by = 0.5),1))  

  pAB<-p+ annotate("text", x = 4, y = 1, label = "A")+
    annotate("text", x = 1, y = 4, label = "B")
  
  return(list(p=p,
              pAB=pAB))

  
}

metastudies_plot<-function(X,sigma,sign){
  
  significant<-(abs(X/sigma)>critval);
  nooutlier<-sigma<50;
  dat<-data.frame(xvar=X,
                  yvar=sigma,
                  significant = as.factor(significant&nooutlier))
  rangeX=1.1*max(max(abs(X)), max(abs(sigma[nooutlier]))*critval);
  
  dat<-dat[order(dat$significant),]
  
  if (sign == -1 ) {
    p<-ggplot(dat, aes(x=xvar,y=yvar)) +
      xlim(c(-rangeX,rangeX))+
      ylim(c(0,rangeX/critval))+
      xlab("X")+
      ylab(expression("sigma" ))+
      geom_abline(intercept = 0,slope=1/critval,color="grey")+ 
      geom_abline(intercept = 0,slope=-1/critval,color="grey")+
      geom_point(shape=shape,size = 2,aes(colour = significant,
                                          fill = significant))+ 
      scale_fill_manual(values=c("grey", "blue")) + 
      scale_colour_manual(values=c("grey", "blue"))+
      scale_y_continuous( expand = c(0,0), limits = c(0,rangeX/critval))
# scale_y_continuous( expand = c(0,0), limits = c(0,rangeX/critval),
#                        breaks = round(seq(min(dat$yvar), max(dat$yvar), by = 0.5),1))+
#      scale_x_continuous(breaks = round(seq(min(dat$xvar), max(dat$xvar), by = 0.5),1))  
    
  } else {
    p<-ggplot(dat, aes(x=xvar,y=yvar)) +
      xlim(c(0,rangeX))+
      ylim(c(0,rangeX/critval))+
      xlab("|X|")+
      ylab(expression("sigma" ))+
      geom_abline(intercept = 0,slope=1/critval,color="grey")+ 
      geom_point(shape=shape,size = 2,aes(colour = significant,
                                          fill = significant))+ 
      scale_fill_manual(values=c("grey", "blue")) + 
      scale_colour_manual(values=c("grey", "blue"))+
      scale_y_continuous( expand = c(0,0), limits = c(0,rangeX/critval))+
      scale_x_continuous( expand = c(0,0), limits = c(0,rangeX))
#    ,
#                          breaks = round(seq(min(dat$yvar), max(dat$yvar), by = 0.5),1))+
#      scale_x_continuous(breaks = round(seq(min(dat$xvar), max(dat$xvar), by = 0.5),1))
  }
  
 
  
  return(p)
}
do_hist<-function(Zuse,symmetric,xaxisname) {

  ll=floor(min(Zuse));
  lleven=floor(ll/2)*2;
  uu=ceiling(max(Zuse));
  
  if (symmetric == 0) {
    if (number>=30) {
      uu2<-ceiling((uu-.36)/.32)*.32+.36;
      ll2<-floor((ll+.36)/.32)*.32-.36;
      edges<-c(seq(from=ll2,
                   to=-0.36,
                   by=0.32), 0, seq(from=0.36,
                                    to=uu2,
                                    by=0.32));
    } else {
      uu2<-ceiling((uu-.68)/.64)*.64+.68;
      ll2<-floor((ll+.68)/.64)*.64-.68;
      edges<-c(seq(from=ll2,
                   to=-0.68,
                   by=0.64), 0, seq(from=0.68,
                                    to=uu2,
                                    by=0.64));
    }
  } else if (symmetric==1){
    if (number>=30) {
      uu2=ceiling((uu-.36)/.32)*.32+.36;
      edges = c(0,seq(from=.36,
                      to=1.64,
                      by=0.32),
                seq(from=1.96,
                    to=uu2,
                    by=0.32)
      )
    }  else {
      uu2=ceiling((uu-.68)/.64)*.64+.68;
      edges= c(0, seq(from=.68,
                      to=uu2,
                      by=0.64))
    }
    
  }
  
  
  
  if (symmetric == 0) {
    
    h<-ggplot(data = as.data.frame(Zuse), aes(Zuse))+
      geom_histogram(aes(y = ..density..),
                     fill = 'blue',
                     breaks=edges)+
      geom_vline(xintercept =-1.96,color='grey')+
      geom_vline(xintercept =1.96, color='grey')+
      xlab(xaxisname)+
      ylab('Density')+
      xlim(c(min(edges),max(edges)))
#      scale_x_continuous(breaks = round(seq(min(Zuse), max(Zuse), by = 0.5),1))
    
    
  } else {
    h<-ggplot(data = as.data.frame(Zuse), aes(Zuse))+
      geom_histogram(aes(y = ..density..),
                     fill = 'blue',
                     breaks=edges)+
      geom_vline(xintercept =1.96, color='grey')+
      xlab(xaxisname)+
      ylab('Density')+
      xlim(c(min(edges),max(edges)))
#      scale_x_continuous(breaks = round(seq(min(Zuse), max(Zuse), by = 0.5),1))
  }
  return(h)
}

