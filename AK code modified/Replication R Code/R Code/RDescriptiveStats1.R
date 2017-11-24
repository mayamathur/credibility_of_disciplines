RDescriptiveStats1<-function(X,X_meta,sigma_meta,combined,name,symmetric)  {

  show(application)
  source("Plot_functions.R")
  
  if (combined) {
    X_rep<-X
    
    xaxisname<-"W"
    yaxisname<-"WRep"
    A<-replication_plot(X_rep,xaxisname,yaxisname)
    p<-A$p
    pAB<-A$pAB
    
   
    
    p_meta<-metastudies_plot(X_meta,sigma_meta,-1)
    Zuse<-X_meta/sigma_meta;
    show(xaxisname)
    p_hist<-do_hist(Zuse,symmetric,xaxisname);
    
    
    p12<-grid.arrange(p, p_meta, ncol = 2)
    filepath<-paste0(pathname,'/FiguresandTables/',name, 'CombinedScatter.pdf')
    save_plot(filepath,p12,ncol=2,base_width = 4, base_height=3)
    dev.off()
    
    p_meta<-metastudies_plot(abs(X_meta),sigma_meta,1)
    p123<-grid.arrange(p_hist,p, p_meta, ncol = 3)
    filepath<-paste0(pathname,'/FiguresandTables/',name, 'CombinedScatterHist.pdf')
    save_plot(filepath,p123,ncol=3,base_width = 4, base_height=3)
    dev.off()
    show(name)
  } else {
    xaxisname<-"Z"
    yaxisname<-"ZRep"
    
    #show(xaxisname)
    
    A<-replication_plot(X,xaxisname,yaxisname)
    p<-A$p
    pAB<-A$pAB
    
    filepath<-paste0(pathname,'/FiguresandTables/',name, 'Scatter.pdf')
    pdf(filepath, width=5, height=5)
    print(p)
    dev.off()
    
    filepath<-paste0(pathname,'/FiguresandTables/',name, 'SlidesScatter.pdf')
    pdf(filepath, width=5, height=5)
    print(pAB)
    dev.off()
    
  }
  
}