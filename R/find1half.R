#' find.1half()
#'
#' Analytic algorithm for finding the OPEE half-life and effect size.
#' @param h0 Initialized half-life parameter to use in the OPEE estimation process. Default is NULL, though an error will be returned if neither h0 nor k0 are provided.
#' @param k0 Initial rate to use for effective exposure estimation process. Default is NULL. Preferable to input h0 instead, though algorithm uses k0.
#' @param datter Dataframe with the corresponding columns needed in the estimation algorithm process. Please see additional documentation for appropriate structure of the dataset, including the required column names and types.
#' @param modtype Whether to use survival or pooled logistic regression methods. Default is \code{modtype="cph"} indicating to use survival. If anything else, will assume "pool".
#' @param stdose Specify the default standard dose for which the effective exposure plateau is desired. Default assumes that this is a binary exposure with value of 1 implying exposure periods.
#' @param covar Whether or not to adjust for covariates during the half-life findind process. Input expects that these will match the variable names (columns of the dataframe object) that are to be use. Default is NULL.
#' @param stratas For survival outcomes that require separate strata, indicate the list of stratas here. Input expects that these will match the variable names (columns of the dataframe object). Default is NULL.
#' @param iterh Indicator for whether or not to track all iterations as the algorithm progresses
#' @param printer Indicator for whether or not to output (in the console) the iterations of the algorithm
#' @param txtadd Text to add onto print-out (if specified by \code{printer=1}). Used for simulations tracking during parallel processing.
#' 
#' @return Multiple objects related to the algorithmic convergence and final result using a one-parameter effective exposure model. (More details to be added)
#' @export
find.1half=function(h0=NULL,k0=NULL,datter,modtype="cph",stdose=1,covar=NULL,stratas=NULL,iterh=NULL,printer=NULL,txtadd=NULL){
  
  ####
  ###
  #       data preparedness checks
  ###
  ####
  
  if(is.null(k0) & is.null(h0)) return(print("ERROR: Specify either k0 or h0"))
  
  if(!is.null(k0) & is.null(h0)){
    h.0<-log(2)/k0
  } else {
    h.0=h0
  }
  
  # Specify the column pointers for the most-used data vectors
  t1<-grep("time",names(datter),fixed=TRUE)
  t2<-grep("tstop",names(datter),fixed=TRUE)
  case<-grep("event",names(datter),fixed=TRUE)
  
  if(length(covar)>0) {
    cv<-which(names(datter) %in% covar)
  } else {
    cv<-NULL
  }
  
  if(length(stratas)>0) {
    s<-which(names(datter) %in% stratas)
  } else {
    s<-NULL
  }
  
  ####
  ###
  #       initial model fits
  ###
  ####
  
  hgrid<-c(h.0/2,h.0,h.0*2)  # initialize the first three half-lives
  
  # Calculate the effective exposure concentration for each half-life
  Cgrid<-sapply(1:3,function(i) C1fun.h(thalf=hgrid[i],dat=datter))
  
  # Fit the half-life specific models
  if(modtype=="cph") {
    mgrid<-lapply(1:3,function(i) quick.cph(C.in=Cgrid[,i],dat=datter,t1=t1,t2=t2,case=case,strats=s,covs=cv))
    # Pull out the model components - log likelihood, beta1, Ratio
    llgrid<-sapply(1:3,function(i) mgrid[[i]]$logLik)
    b1.grid<-sapply(1:3, function(c) mgrid[[c]]$beta.C)
  } else {
    mgrid<-lapply(1:3,function(i) fit.pool(C.in=Cgrid[,i],dat=datter,covs=covar))
    # Pull out the model components - log likelihood, beta1, Ratio
    llgrid<-sapply(1:3,function(i) as.numeric(logLik(mgrid[[i]])))
    b1.grid<-sapply(1:3, function(c) as.numeric(mgrid[[c]]$coefficients["C.in"]))
  }
  
  Ratiogrid<-sapply(1:3, function(q) exp(b1.grid[q]*stdose))
  
  rm(Cgrid,mgrid)
  
  # Set the initial difference markers to NA
  diff.h<-NA
  diff.Ratio<-NA
  diff.llike<-NA
  algostep<-0
  numiter<-0
  
  half<-hgrid[2]
  llike<-llgrid[2]
  beta1<-b1.grid[2]
  Ratioest<-Ratiogrid[2]
  
  if(length(printer)>0) {
    init.fit<-c(round(half,2),round(beta1,3),round(llike,2))
    addtxt<-ifelse(length(txtadd)>0,paste("Sim:",paste(txtadd,c(", Starting Values (Lag Half-Life, Beta, LogLik)"))),c("Starting Values (Lag Half-Life, Beta, LogLik)"))     
    print(addtxt)
    print(init.fit)
  }
  
  if(length(iterh)>0){
    # Output the initialized half-life values
    iter.hist<-data.frame(diff.h=diff.h,
                          half=half,
                          llike=round(llike,3),
                          diff.llike=diff.llike,
                          beta1=round(beta1,3),
                          Ratioest=round(Ratioest,3),
                          diff.Ratio=diff.Ratio,
                          step=algostep)
  }
  
  # Combine the logL results with the half-life params
  iterpath<-as.data.frame(cbind(hgrid,llgrid,b1.grid))
  iterpath$id<-rownames(iterpath)
  
  if(length(printer)>0) {
    addtxt<-ifelse(length(txtadd)>0,paste("Sim:",paste(txtadd,c(", Step:"))),c("Step:"))     
    print(paste(addtxt,numiter))
    print(algostep)
    print(iterpath)
  }
  iterout<-data.frame(iterpath,step=numiter)
  
  # Identify the maximum logL for the first 3 combinations
  index0<-as.numeric(iterpath$id[which.max(iterpath$llgrid)])
  if(length(index0)>0){
    index<-ifelse(any(index0==2),2,ifelse(all(index0>=2),min(index0),max(index0)))
    rm(index0)
  }else {index<-index0}
  
  
  
  ####
  ###
  #       First iterations - 
  #         Drop or increase to the area of the max - i.e. center-up
  ###
  ####
  
  while(index != 2 & numiter < 100 & (diff.llike>1e-2 | is.na(diff.llike))){
    algostep<-1
    numiter<-numiter+1
    
    llike2<-llgrid[index]
    half2<-hgrid[index]
    Ratioest2<-Ratiogrid[index]
    beta2<-b1.grid[index]
    
    diff.h<-abs(half2-half) # determine the distance between index h and center
    diff.llike<-abs(llike2-llike) # determine the distance between max logLik and center
    diff.Ratio<-abs(Ratioest2-Ratioest)  # determine the distance between index Ratio and center
    
    if(length(iterh)>0){
      # append the summary history
      iter.hist<-rbind(iter.hist,
                       cbind(diff.h=round(diff.h,2),
                             half=half2,
                             llike=round(llike2,3),
                             diff.llike=round(diff.llike,3),
                             beta1=round(beta2,3),
                             Ratioest=round(Ratioest2,3),
                             diff.Ratio=round(diff.Ratio,3),
                             step=algostep))
    }
    
    if(index==1){
      hgrid[2:3]<-hgrid[1:2]
      # Cgrid[,2:3]<-Cgrid[,1:2]
      llgrid[2:3]<-llgrid[1:2]
      b1.grid[2:3]<-b1.grid[1:2]
      Ratiogrid[2:3]<-Ratiogrid[1:2]
      
      hgrid[1]<-hgrid[2]/2
      Cnew<-C1fun.h(thalf=hgrid[1],dat=datter)
    } else if(index==3){
      hgrid[1:2]<-hgrid[2:3]
      # Cgrid[,1:2]<-Cgrid[,2:3]
      llgrid[1:2]<-llgrid[2:3]
      b1.grid[1:2]<-b1.grid[2:3]
      Ratiogrid[1:2]<-Ratiogrid[2:3]
      
      hgrid[3]<-hgrid[2]*2
      Cnew<-C1fun.h(thalf=hgrid[3],dat=datter)
    }
    
    if(modtype=="cph") {
      newmod<-quick.cph(C.in=Cnew,dat=datter,t1=t1,t2=t2,case=case,strats=s,covs=cv)
      llgrid[index]<-newmod$logLik
      b1.grid[index]<-newmod$beta.C
    } else {
      newmod<-fit.pool(C.in=Cnew,dat=datter,covs=covar)
      llgrid[index]<-as.numeric(logLik(newmod))
      b1.grid[index]<-as.numeric(newmod$coefficients["C.in"])
    }
    
    rm(Cnew, newmod)
    Ratiogrid[index]<-exp(b1.grid[index]*stdose)
    
    iterpath<-as.data.frame(cbind(hgrid,llgrid,b1.grid))
    iterpath<-iterpath[order(iterpath$hgrid),]
    iterpath$id<-rownames(iterpath)
    if(length(printer)>0) {
      addtxt<-ifelse(length(txtadd)>0,paste("Sim:",paste(txtadd,c(", Step:"))),c("Step:"))          
      print(paste(addtxt,numiter)) 
      print(algostep)
      print(iterpath)
    }
    iterout<-rbind(iterout,data.frame(iterpath,step=numiter))
    llike<-llike2
    beta1<-beta2
    Ratioest<-Ratioest2
    half<-half2
    
    index0<-as.numeric(iterpath$id[which.max(iterpath$llgrid)])
    if(length(index0)>1){
      index<-ifelse(any(index0==2),2,ifelse(all(index0>=2),min(index0),max(index0)))
      rm(index0)
    }else {index<-index0}
    
  }
  
  ####
  ###
  #  Second iterations - start once centered value is max loglik
  # Adjust the right side to the same distance from center as left
  # At this point, the algorithm sets [1] as 1/2 of [2] which is 1/2 of [3]
  #         Reset [3] to be 3*[1]
  #         Iterate through halving distance until center is no longer [2]
  ###
  ####
  
  algostep<-2
  numiter<-numiter+1
  # Set the distance marker
  diff.h<-hgrid[1]
  
  # Calculate the new half [3] info
  hgrid[3]<-3*hgrid[1]
  
  C.3<-C1fun.h(thalf=hgrid[3],dat=datter)
  
  if(modtype=="cph") {
    m.3<-quick.cph(C.in=C.3,dat=datter,t1=t1,t2=t2,case=case,strats=s,covs=cv)
    llgrid[3]<-m.3$logLik
    b1.grid[3]<-m.3$beta.C
  } else {
    m.3<-fit.pool(C.in=C.3,dat=datter,covs=covar)
    llgrid[3]<-as.numeric(logLik(m.3))
    b1.grid[3]<-as.numeric(m.3$coefficients["C.in"])
  }
  
  rm(C.3, m.3)
  Ratiogrid[3]<-exp(b1.grid[3]*stdose)
  iterpath<-as.data.frame(cbind(hgrid,llgrid,b1.grid))
  iterpath<-iterpath[order(iterpath$hgrid),]
  iterpath$id<-rownames(iterpath)
  if(length(printer)>0) {
    addtxt<-ifelse(length(txtadd)>0,paste("Sim:",paste(txtadd,c(", Step:"))),c("Step:"))          
    print(paste(addtxt,numiter)) 
    print(algostep)
    print(iterpath)
  }
  iterout<-rbind(iterout,data.frame(iterpath,step=numiter))
  
  index0<-as.numeric(iterpath$id[which.max(iterpath$llgrid)])
  if(length(index0)>1){
    index<-ifelse(any(index0==2),2,ifelse(all(index0>=2),min(index0),max(index0)))
    rm(index0)
  }else {index<-index0}
  
  while(index==2 & (diff.llike > 1e-2 | is.na(diff.llike)) & numiter < 100){
    numiter<-numiter+1
    ll2<-max(llgrid[-2]) # Find the second max logLike
    i2<-which(llgrid==ll2)
    
    diff.llike<-abs(llgrid[i2]-llgrid[2])
    diff.Ratio<-abs(Ratiogrid[i2]-Ratiogrid[2])
    
    if(length(iterh)>0){
      iter.hist <- rbind(iter.hist, 
                         data.frame(diff.h=round(diff.h,2), 
                                    half=hgrid[2], 
                                    llike=round(llgrid[2],3), 
                                    diff.llike=round(diff.llike,3),
                                    beta1=round(b1.grid[2],3),
                                    Ratioest=round(Ratiogrid[2],3),
                                    diff.Ratio=round(diff.Ratio,3),
                                    step=algostep))    
    }
    
    # Set the distance marker
    diff.h<-diff.h/2
    hgrid[1]<-hgrid[2]-diff.h
    hgrid[3]<-hgrid[2]+diff.h
    
    C.1<-C1fun.h(thalf=hgrid[1],dat=datter)
    C.3<-C1fun.h(thalf=hgrid[3],dat=datter)
    
    if(modtype=="cph") {
      m.1<-quick.cph(C.in=C.1,dat=datter,t1=t1,t2=t2,case=case,strats=s,covs=cv)
      llgrid[1]<-m.1$logLik
      b1.grid[1]<-m.1$beta.C
      m.3<-quick.cph(C.in=C.3,dat=datter,t1=t1,t2=t2,case=case,strats=s,covs=cv)
      llgrid[3]<-m.3$logLik
      b1.grid[3]<-m.3$beta.C
    } else {
      m.1<-fit.pool(C.in=C.1,dat=datter,covs=covar)
      llgrid[1]<-as.numeric(logLik(m.1))
      b1.grid[1]<-as.numeric(m.1$coefficients["C.in"])
      m.3<-fit.pool(C.in=C.3,dat=datter,covs=covar)
      llgrid[3]<-as.numeric(logLik(m.3))
      b1.grid[3]<-as.numeric(m.3$coefficients["C.in"])
    }
    
    rm(C.3, C.1, m.3, m.1)
    Ratiogrid[1]<-exp(b1.grid[1]*stdose)
    Ratiogrid[3]<-exp(b1.grid[3]*stdose)
    
    iterpath<-as.data.frame(cbind(hgrid,llgrid,b1.grid))
    iterpath<-iterpath[order(iterpath$hgrid),]
    iterpath$id<-rownames(iterpath)
    if(length(printer)>0) {
      addtxt<-ifelse(length(txtadd)>0,paste("Sim:",paste(txtadd,c(", Step:"))),c("Step:"))          
      print(paste(addtxt,numiter)) 
      print(algostep)
      print(iterpath)
    }
    iterout<-rbind(iterout,data.frame(iterpath,step=numiter))
    
    index0<-as.numeric(iterpath$id[which.max(iterpath$llgrid)])
    if(length(index0)>1){
      index<-ifelse(any(index0==2),2,ifelse(all(index0>=2),min(index0),max(index0)))
      rm(index0)
    } else {index<-index0}
    
  }
  
  diff.llike<-abs(llgrid[index]-llgrid[2])
  diff.Ratio<-abs(Ratiogrid[index]-Ratiogrid[2])
  
  if(length(iterh)>0){
    iter.hist <- rbind(iter.hist, 
                       cbind(diff.h=round(diff.h,2), 
                             half=hgrid[index],
                             llike=round(llgrid[index],3), 
                             diff.llike=round(diff.llike,3),
                             beta1=round(b1.grid[index],3), 
                             Ratioest=round(Ratiogrid[index],3),
                             diff.Ratio=round(diff.Ratio,3),
                             step=algostep))   
  }
  
  ####
  ###
  #  Third iterations - Center direction is determined
  #  Adjust the pointer between the two maxima until sufficient stabilization in loglik
  ###
  ####
  
  while(diff.llike > 1e-2 & numiter < 100){
    numiter<-numiter+1
    algostep<-3
    
    ll2<-max(llgrid[-index]) # Find the second max logLike
    i2<-which(llgrid==ll2)
    
    h.1<-mean(c(hgrid[index],hgrid[i2])) # set a new half-life to half-way between the two max points
    diff.h<-abs(h.1-hgrid[index])
    
    C.1<-C1fun.h(thalf=h.1,dat=datter)  
    
    if(modtype=="cph") {
      m.1<-quick.cph(C.in=C.1,dat=datter,t1=t1,t2=t2,case=case,strats=s,covs=cv)
      llike.1<-m.1$logLik
      b1.1<-m.1$beta.C
    } else {
      m.1<-fit.pool(C.in=C.1,dat=datter,covs=covar)
      llike.1<-as.numeric(logLik(m.1))
      b1.1<-as.numeric(m.1$coefficients["C.in"])
    }
    
    rm(C.1, m.1)
    
    Ratio.1<-exp(b1.1*stdose)
    diff.llike<-abs(llike.1-llgrid[index])
    diff.Ratio<-abs(Ratio.1-Ratiogrid[index])
    
    hgrid<-c(hgrid[index],h.1,hgrid[i2])
    llgrid<-c(llgrid[index],llike.1,llgrid[i2]); llgrid<-llgrid[order(hgrid)]
    b1.grid<-c(b1.grid[index],b1.1,b1.grid[i2]); b1.grid<-b1.grid[order(hgrid)]
    Ratiogrid<-c(Ratiogrid[index],Ratio.1,Ratiogrid[i2]); Ratiogrid<-Ratiogrid[order(hgrid)]
    hgrid<-hgrid[order(hgrid)]
    
    iterpath<-as.data.frame(cbind(hgrid,llgrid,b1.grid))
    iterpath<-iterpath[order(iterpath$hgrid),]
    rownames(iterpath)<-NULL
    iterpath$id<-rownames(iterpath)
    if(length(printer)>0) {
      addtxt<-ifelse(length(txtadd)>0,paste("Sim:",paste(txtadd,c(", Step:"))),c("Step:"))          
      print(paste(addtxt,numiter)) 
      print(algostep)
      print(iterpath)
    }
    iterout<-rbind(iterout,data.frame(iterpath,step=numiter))
    
    if(length(iterh)>0){
      iter.hist <- rbind(iter.hist, 
                         cbind(diff.h=round(diff.h,2), 
                               half=h.1, 
                               llike=round(llike.1,3), 
                               diff.llike=round(diff.llike,3),
                               beta1=round(b1.1,3),
                               Ratioest=round(Ratio.1,3), 
                               diff.Ratio=round(diff.Ratio,3),
                               step=algostep))    
    }
    
    index0<-as.numeric(iterpath$id[which.max(iterpath$llgrid)])
    if(length(index0)>1){
      index<-ifelse(any(index0==2),2,ifelse(all(index0>=2),min(index0),max(index0)))
      rm(index0)
    } else {index<-index0}
  }
  
  h<-hgrid[index]
  Cfin<-C1fun.h(thalf=h,dat=datter)
  
  # Refit final model for CoxPH with complete info
  if(modtype=="cph"){
    if(!is.null(stratas)) {
      stratnew<-sapply(1:length(stratas),function(t) paste0("strata(",paste0(stratas[t],")")))
      covar2<-c(covar,stratnew)
    } else {
      covar2<-covar
    }
    m.1<-fit.cph(C.in=Cfin,dat=datter,covs=covar2)
    llike.1<-as.numeric(logLik(m.1))
    b1.1<-as.numeric(m.1$coefficients["C.in"])
  } else {
    m.1<-fit.pool(C.in=Cfin,dat=datter,covs=covar)
    llike.1<-as.numeric(logLik(m.1))
    b1.1<-as.numeric(m.1$coefficients["C.in"])
  }
  
  llike<-round(llike.1,4)
  AIC<-2*(length(m.1$coefficients)+1)-2*llike # AIC needs additional parameter for half-life estimation
  b1<-b1.1
  
  b1.se<-ifelse(modtype=="cph",
                summary(m.1)$coefficients["C.in","se(coef)"],
                summary(m.1)$coefficients["C.in","Std. Error"])
  
  Ratio<-round(exp(b1*stdose),4)
  Ratio.95L<-round(exp(b1*stdose-1.96*b1.se*stdose),4)
  Ratio.95U<-round(exp(b1*stdose+1.96*b1.se*stdose),4)
  
  out <- list()
  if(length(iterh)>0){ 
    out$iter.hist<-iter.hist
    out$niter<-dim(iter.hist)[1]
  }
  
  out$final.results<-c(Beta1.Est=round(b1,4),
                       Beta1.SE=round(b1.se,4),
                       Ratio=Ratio,
                       Ratio.95L=Ratio.95L,
                       Ratio.95U=Ratio.95U,
                       Half.in=round(h,2), 
                       Half.out=round(h,2), 
                       h0=h.0,
                       logL=llike,
                       AIC=round(AIC,4),
                       numiter=numiter,
                       modtype=modtype,
                       estalgo="1param",
                       RatioDose=stdose)
  out$final.mod<-summary(m.1)
  out$extra<-iterout
  
  return(out)
}