source("Plot_functions.R")
RDescriptiveStats2<-function(X, sigma, name,symmetric,cluster_ID) {
  
    critval<-1.96;
    p<-metastudies_plot(X,sigma,-1)
    filepath<-paste0(pathname,'/FiguresandTables/',name, 'Scatter.pdf')
    pdf(filepath, width=5, height=5)
    print(p)
    dev.off()
    
    
    filepath<-paste0(pathname,'/FiguresandTables/',name, 'SlidesScatter.pdf')
    pdf(filepath, width=5, height=5)
    print(p)
    dev.off()
    
    Zuse = X/sigma;
    if (max(abs(Zuse))>10) {
      Zuse<-Zuse[abs(Zuse)<6];
    } 
    h<-do_hist(Zuse,symmetric,'X/sigma')
    # handp<-multiplot(h,p,cols=2)
    handp<-grid.arrange(h, p, ncol = 2)
    filepath<-paste0(pathname,'/FiguresandTables/',name, 'ScatterHist.pdf')
    save_plot(filepath,handp,ncol=2,base_width = 4, base_height=3)
  }
  
  
  #  %Meta-regression estimates
  if (identificationapproach==2 && symmetric==0) {
    source("MetaRegressionTable.R")
    source("Clustered_covariance_estimate.R")
    do_metatable<-function(sigma,Zstats) {
      R = cbind(rep(1,length(X)),sigma);
      betahat<-solve(t(R)%*%R)%*%t(R)%*%X
      ehat<-X-R%*%betahat
      
      Sigma_base<-R*matrix(rep(ehat,dim(R)[2]),ncol=dim(R)[2],byrow=FALSE)
      Sigma<-Clustered_covariance_estimate(Sigma_base,cluster_ID);
      Vhat<-number*solve((t(R)%*%R))%*%Sigma%*%solve((t(R)%*%R))
      se_robust<-sqrt(diag(Vhat))
      
      betahat<-as.vector(betahat)
      MetaRegressionTable(pathname,betahat,se_robust,name,Zstats)
    }
    
    do_metatable(sigma,0)
    
    do_metatable(sigma^(-1),1)
    
  
  
}