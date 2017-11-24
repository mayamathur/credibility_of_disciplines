output_result<-function(res,
                        true.value,filename,directoryname,
                        methods,add_info,
                        alpha = 0.05,
                        digs=3, measures=c("Bias","St.Error","RMSE","Rej.Freq.","MedianBias","Range"),
                        N_rep) {
  
  b.hat<-res$b.hat
  st.error.hat<-res$st.error.hat
  norm.b.hat<-res$norm.b.hat
  which_coefs<-res$which_coefs
  #show(b.hat)
  
  Results<-array(0,c(length(true.value),length(methods),length(measures))) 
  dimnames(Results)<-list(names(true.value),methods,measures)
  for (n in methods) {
    #est.name<-strsplit(n,split=split.symb)[[1]][1]
    #mode.name<-strsplit(n,split=split.symb)[[1]][2]
    show(n)
    b<-as.matrix(b.hat[[n]][which_coefs,1:N_rep])
    # show(b)
    st.error<-st.error.hat[[n]][which_coefs,1:N_rep]
    rej.freq<-as.matrix(norm.b.hat[[n]][which_coefs,1:N_rep])
    Results[,n,"Bias"]<-apply(b,1,mean)-true.value
    Results[,n,"St.Error"]<-apply(b,1,sd)
    Results[,n,"RMSE"]<-sqrt((Results[,n,"Bias"])^2+(Results[,n,"St.Error"])^2)
    Results[,n,"Rej.Freq."]<-as.numeric(apply(abs(rej.freq)>qnorm(1-alpha/2),1,mean))
    #show(Results[,n,"Rej.Freq."])
    #show(n)
    Results[,n,"MedianBias"]<-apply(b,1,median)-true.value
    Results[,n,"Range"]<-(apply(b,1,quantile, probs=0.95)- apply(b,1,quantile, probs=0.05))/3.28
    
  #  Results[,n,"Av.Conf.Int.Length"]<-as.numeric(apply(2*st.error,1,mean))
  }
  add_info_caption<-c()
  for (j in 1:length(add_info)) {
    add_info_caption<-c(add_info_caption,paste (names(add_info)[j],add_info[j][[1]]))
    
  }
  add_info_caption<-paste0(add_info_caption,collapse=".")  
  
  
  Results<-t(apply(Results,1,as.vector))
  rownames(Results)<-dimnames(Results)[[1]]
  colnames(Results)<-rep(methods,length(measures))
  caption<-paste0(paste0(measures,collapse=", "),'.',
                  paste0(names(methods),collapse=", "), '.', 
                  add_info_caption)
  label<-paste0(paste0(measures,collapse=", "))
  #filename<-"/Users/virasemenora/Dropbox (MIT)/18.657/Tex output/LROLS.tex"
  setwd("..")
  setwd(directoryname)
  to.Latex(digs=digs,Results,cap=caption,lab=label,filename=filename,k=length(methods),step=length(measures),measures=measures)
  
  
}

#### Auxiliary functions
#### Create latex table "Bias,St.error,RMSE,Rej.freq. of different estimators
to.Latex<-function (digs,matrix,cap,lab,filename,k,step=4,measures=c("Bias","St.Error","RMSE","Rej.Freq.")) {
  #options(digits=12)
  align<-rep("r",k)
  align<-paste0(align,collapse="")
  align<-paste0(align,"|")
  align<-rep(align,step)
  align<-paste0(align,collapse="")
  align<-paste0("r|",align)
  tab<-xtable(matrix,digits=digs,caption=cap,label=lab,align=align)
  #caption(tab)<-cap
  #label(tab)<-lab
  addtorow <- list()
  addtorow$pos <- list(0)
  addtorow$command <- paste0(paste0('& \\multicolumn{',as.character(k),'}{c}{', measures, '}', collapse=''), '\\\\')
  show(digs)
  print(tab,file=filename,append=TRUE, 
        include.rownames=TRUE,
        include.colnames=TRUE,
        sanitize.text.function=function(x){x},
        add.to.row=addtorow)
  
  
}
## Raise vector vec to power p
pow<-function(vec,p) {
  return(vec^p)
}
## Second norm of norm_vec
norm_vec <- function(x) {sqrt(sum(x^2))}

# Compute Robust ST ERRORS
robust_se<-function(Psihat,stepsize) {
  Info = matrix(0,length(Psihat),length(Psihat));
  for (n1 in 1:length(Psihat)) {
    for (n2 in 1:length(Psihat)) {
      thetaplusplus=Psihat;
      thetaplusminus=Psihat;
      thetaminusplus=Psihat;
      thetaminusminus=Psihat;
      
      thetaplusplus[n1]=thetaplusplus[n1]+stepsize;
      thetaplusplus[n2]=thetaplusplus[n2]+stepsize;
      LLH_plusplus=LLH_only(thetaplusplus);
      
      thetaplusminus[n1]=thetaplusminus[n1]+stepsize;
      thetaplusminus[n2]=thetaplusminus[n2]-stepsize;
      LLH_plusminus=LLH_only(thetaplusminus);
      
      thetaminusplus[n1]=thetaminusplus[n1]-stepsize;
      thetaminusplus[n2]=thetaminusplus[n2]+stepsize;
      LLH_minusplus=LLH_only(thetaminusplus);
      
      thetaminusminus[n1]=thetaminusminus[n1]-stepsize;
      thetaminusminus[n2]=thetaminusminus[n2]-stepsize;
      LLH_minusminus=LLH_only(thetaminusminus);
      
      
      Info[n1,n2]=((LLH_plusplus-LLH_plusminus)/(2*stepsize)-(LLH_minusplus-LLH_minusminus)/(2*stepsize))/(2*stepsize);
      
    }
  }
  Var = solve(Info)
  se = sqrt(diag(Var))
  return(se)
}