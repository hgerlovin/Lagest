#' find.2half()
#'
#' Analytic algorithm for finding the TPEE half-life and effect size.
#' @inheritParams find.1half
#' @param h0 Same as with find.1half, however, the assumed parameter pair is set to the single half-life at the start of the TPEE search algorithm.
#' @param tol.in Granularity of convergence at which to stop looking for the incline half-life parameter. Default is 1e-2.
#' @param tol.out Granularity of convergence at which to stop looking for the decline half-life parameter. Default is 1e-2.
#' @param tol.ll At which point iteration towards maximum negative log-likelihood can stop. Default is 1e-3.
#' 
#' @return Multiple objects related to the algorithmic convergence and final result using a two-parameter effective exposure model. (More details to be added)
#' @export
find.2half=function(h0=NULL,k0=NULL,modtype="cph",datter,stdose=1,covar=NULL,stratas=NULL,tol.in=1e-2,tol.out=1e-2,tol.ll=1e-3,iterh=NULL,printer=NULL,txtadd=NULL){
  
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
  
  hgrid<-c(h.0/2,h.0,h.0*2)  # initialize the first three half-lives
  
  h.in<-hgrid
  h.out<-hgrid
  
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
  
  # Create a grid of all 9 initialized values
  hgrid0<-expand.grid(h.in=h.in,h.out=h.out)
  
  # Calculate the corresponding concentrations
  Cgrid0<-sapply(1:9,function(i) C1fun.2h(thalf=c(hgrid0$h.in[i],hgrid0$h.out[i]),dat=datter)$Conc)
  
  # Fit the initial models
  if(modtype=="cph") {
    mgrid0<-lapply(1:9,function(i) quick.cph(C.in=Cgrid0[,i],dat=datter,t1=t1,t2=t2,case=case,strats=s,covs=cv))
    # Pull out the model components - log likelihood, beta1, Ratio
    llgrid0<-sapply(1:9,function(i) mgrid0[[i]]$logLik)
    beta1.0<-sapply(1:9, function(c) mgrid0[[c]]$beta.C)
  } else {
    mgrid0<-lapply(1:9,function(i) fit.pool(C.in=Cgrid0[,i],dat=datter,covs=covar))
    # Save the initial logL's and betas for EE
    llgrid0<-sapply(1:9,function(i) as.numeric(logLik(mgrid0[[i]])))
    beta1.0<-sapply(1:9, function(c) as.numeric(mgrid0[[c]]$coefficients["C.in"]))
  }
  
  rm(Cgrid0,mgrid0)
  numiter<-0
  # Set the initial difference markers to NA
  diff.inhf<-NA
  diff.outhf<-NA
  diff.llike<-NA
  ll1<-llgrid0[5]
  hin1<-hgrid0$h.in[5]
  hout1<-hgrid0$h.in[5]
  beta1<-beta1.0[5] 
  
  if(length(printer)>0) {
    init.fit<-c(round(hin1,2),round(hout1,2),round(beta1,3),round(ll1,2))
    addtxt<-ifelse(length(txtadd)>0,paste("Sim:",paste(txtadd,c(", Starting Values (Incline, Decline, Beta, LogLik)"))),c("Starting Values (Incline, Decline, Beta, LogLik)"))     
    print(addtxt)
    print(init.fit)
  }
  
  if(length(iterh)>0){
    # Output the initialized half-life values
    iter.hist<-data.frame(diff.inhf=diff.inhf,
                          diff.outhf=diff.outhf,
                          half.in=round(hin1,2),
                          half.out=round(hout1,2),
                          llike=round(ll1,2),
                          diff.llike=diff.llike,
                          beta1=round(beta1,3),
                          Ratio=round(exp(beta1*stdose),3),
                          center=paste(round(hin1,2),round(hout1,2),sep="|"))
  }
  
  # Combine the logL results with the half-life params
  iterpath<-cbind(hgrid0,llgrid0,beta1.0)
  iterpath$id<-rownames(iterpath)
  iterpath$hname<-paste(iterpath$h.in,iterpath$h.out,sep="|")
  iterpath<-iterpath[order(iterpath$id),]
  
  if(length(printer)>0) {
    addtxt<-ifelse(length(txtadd)>0,paste("Sim:",paste(txtadd,c(", Step:"))),c("Step:"))          
    print(paste(addtxt,numiter)) 
    print(iterpath)
  }
  iterout<-data.frame(iterpath,step=numiter)
  
  # Identify the maximum logL for the first 9 combinations
  index0<-as.numeric(iterpath$id[which.max(iterpath$llgrid0)])
  if(length(index0)>0){
    index<-ifelse(any(index0==5),5,ifelse(all(index0>=5),min(index0),max(index0)))
    rm(index0)
    index0<-index
    rm(index)
  }
  
  # If the index0 is NOT centered completely, adjust the in/out grid
  
  while(index0 !=5 & numiter < 100 & (diff.llike>tol.ll | is.na(diff.llike))){
    numiter<-numiter+1
    
    # Record the first maximum
    diff.llike<-abs(iterpath$llgrid0[iterpath$id==index0]-ll1)
    diff.inhf<-abs(iterpath$h.in[iterpath$id==index0]-hin1)
    diff.outhf<-abs(iterpath$h.out[iterpath$id==index0]-hout1)
    ll1<-iterpath$llgrid0[iterpath$id==index0]
    beta1<-iterpath$beta1.0[iterpath$id==index0]
    # Update the center values
    hin1<-iterpath$h.in[iterpath$id==index0]
    hout1<-iterpath$h.out[iterpath$id==index0]
    
    if(length(iterh)>0){
      # Output the initialized half-life values
      iter.hist<-rbind(iter.hist,
                       data.frame(diff.inhf=round(diff.inhf,2),
                                  diff.outhf=round(diff.outhf,2),
                                  half.in=round(hin1,2),
                                  half.out=round(hout1,2),
                                  llike=round(ll1,2),
                                  diff.llike=round(diff.llike,3),
                                  beta1=round(beta1,3),
                                  Ratio=round(exp(beta1*stdose),3),
                                  center=paste(round(hin1,2),round(hout1,2),sep="|")))
    }
    
    if(index0==1){
      # index0==1: top left corner, expand to 1/2 of both values
      # shift in-param left
      h.in[2:3]<-h.in[1:2]
      h.in[1]<-h.in[2]/2 
      # shift out-param up 
      h.out[2:3]<-h.out[1:2]
      h.out[1]<-h.out[2]/2 
    } else if(index0==2){
      # index0==2: top middle
      # shift out-param up
      h.out[2:3]<-h.out[1:2]
      h.out[1]<-h.out[2]/2 
      # h.in stays the same
    } else if(index0==3){ 
      # index0==3: top right corner
      # shift in-param right
      h.in[1:2]<-h.in[2:3]
      h.in[3]<-h.in[2]*2 
      # shift out-param up
      h.out[2:3]<-h.out[1:2]
      h.out[1]<-h.out[2]/2 
    } else if(index0==4){ 
      # index0==4: center left side
      # shift in-param left
      h.in[2:3]<-h.in[1:2]
      h.in[1]<-h.in[2]/2 
    } else if(index0==6){ 
      # index0==6: center right side
      # shift in-param right
      h.in[1:2]<-h.in[2:3]
      h.in[3]<-h.in[2]*2 
    } else if(index0==7){ 
      # index0==7: bottom left corner
      # shift in-param left
      h.in[2:3]<-h.in[1:2]
      h.in[1]<-h.in[2]/2 
      # shift out-param down
      h.out[1:2]<-h.out[2:3]
      h.out[3]<-h.out[2]*2 
    } else if(index0==8){ 
      # index0==8: bottom middle
      # shift out-param down
      h.out[1:2]<-h.out[2:3]
      h.out[3]<-h.out[2]*2 
    } else if(index0==9){ 
      # index0==9: bottom right corner
      # shift in-param right
      h.in[1:2]<-h.in[2:3]
      h.in[3]<-h.in[2]*2 
      # shift out-param down
      h.out[1:2]<-h.out[2:3]
      h.out[3]<-h.out[2]*2 
    }
    
    # Setup a new grid with the refined incline and decline parameters
    hnew<-expand.grid(h.in=h.in,h.out=h.out)
    
    # Name the combinations to minimize calculations
    # May want to just pull from the iterpath here instead of simply removing
    hnew$hname<-paste(hnew$h.in,hnew$h.out,sep="|")
    hnew$id<-rownames(hnew)
    
    iterpath<-iterpath[,names(iterpath)!=c("id")]
    tout<-intersect(hnew$hname,iterpath$hname)
    hpull<-merge(hnew,iterpath,by=c("hname","h.in","h.out")) # find way to merge in the logL, h.in, h.out, and beta param only
    hnew<-hnew[!(hnew$hname %in% tout),] # remove the repeated obs
    
    newdim<-length(hnew$h.in)
    
    # Fit the new models
    Cfit<-sapply(1:newdim,function(i) C1fun.2h(thalf=c(hnew$h.in[i],hnew$h.out[i]),dat=datter)$Conc)
    
    if(modtype=="cph") {
      mfit<-lapply(1:newdim,function(i) quick.cph(C.in=Cfit[,i],dat=datter,t1=t1,t2=t2,case=case,strats=s,covs=cv))
      llfit<-sapply(1:newdim,function(i) mfit[[i]]$logLik)
      betafit<-sapply(1:newdim,function(i) mfit[[i]]$beta.C)
    } else {
      mfit<-lapply(1:newdim,function(i) fit.pool(C.in=Cfit[,i],dat=datter,covs=covar))
      llfit<-sapply(1:newdim,function(i) as.numeric(logLik(mfit[[i]])))
      betafit<-sapply(1:newdim,function(i) as.numeric(mfit[[i]]$coefficients["C.in"]))
    }
    
    rm(Cfit,mfit)
    
    newfits<-as.data.frame(cbind(hnew,llgrid0=llfit,beta1.0=betafit))
    
    iterpath<-suppressWarnings(dplyr::bind_rows(hpull,newfits))
    iterpath<-iterpath[order(iterpath$id),]
    rownames(iterpath)<-NULL
    if(length(printer)>0) {
      addtxt<-ifelse(length(txtadd)>0,paste("Sim:",paste(txtadd,c(", Step:"))),c("Step:"))          
      print(paste(addtxt,numiter)) 
      print(iterpath)
    }
    iterout<-rbind(iterout,data.frame(iterpath,step=numiter))
    
    index0<-as.numeric(iterpath$id[which.max(iterpath$llgrid0)])
    if(length(index0)>1){
      index<-ifelse(any(index0==5),5,ifelse(all(index0>=5),min(index0),max(index0)))
      rm(index0)
      index0<-index
      rm(index)
    }
    
  }
  
  numiter<-numiter+1
  ll1<-iterpath$llgrid0[iterpath$id==index0]
  beta1<-iterpath$beta1.0[iterpath$id==index0]
  
  # Now that centered, look at the same distance in the upper directions
  diff.inhf<-iterpath$h.in[iterpath$id==1]
  diff.outhf<-iterpath$h.out[iterpath$id==1]
  
  centin<-iterpath$h.in[iterpath$id==5]
  centout<-iterpath$h.out[iterpath$id==5]
  cent<-paste(round(centin,2),round(centout,2),sep="|")
  
  newin<-centin+diff.inhf
  newout<-centout+diff.outhf
  
  # Create new incline vector that includes two maxes bracketing mean for fixed decline
  hnew.in<-sort(c(diff.inhf,centin,newin))
  # Create new decline vector that includes two maxes bracketing mean for fixed incline
  hnew.out<-sort(c(diff.outhf,centout,newout))
  
  # Setup a new grid with the refined incline and decline parameters
  hnew<-expand.grid(h.in=hnew.in,h.out=hnew.out)
  
  # Name the combinations to minimize calculations
  hnew$hname<-paste(hnew$h.in,hnew$h.out,sep="|")
  hnew$id<-rownames(hnew)
  
  iterpath<-iterpath[,names(iterpath)!=c("id")]
  
  tout<-intersect(hnew$hname,iterpath$hname)
  hpull<-merge(hnew,iterpath,by=c("hname","h.in","h.out")) # find way to merge in the logL, h.in, h.out, and beta param only
  hnew<-hnew[!(hnew$hname %in% tout),] # remove the repeated obs
  
  newdim<-length(hnew$h.in) 
  
  # Fit the new models
  Cfit<-sapply(1:newdim,function(i) C1fun.2h(thalf=c(hnew$h.in[i],hnew$h.out[i]),dat=datter)$Conc)
  
  if(modtype=="cph") {
    mfit<-lapply(1:newdim,function(i) quick.cph(C.in=Cfit[,i],dat=datter,t1=t1,t2=t2,case=case,strats=s,covs=cv))
    llfit<-sapply(1:newdim,function(i) mfit[[i]]$logLik)
    betafit<-sapply(1:newdim,function(i) mfit[[i]]$beta.C)
  } else {
    mfit<-lapply(1:newdim,function(i) fit.pool(C.in=Cfit[,i],dat=datter,covs=covar))
    llfit<-sapply(1:newdim,function(i) as.numeric(logLik(mfit[[i]])))
    betafit<-sapply(1:newdim,function(i) as.numeric(mfit[[i]]$coefficients["C.in"]))
  }
  
  rm(Cfit,mfit)
  
  newfits<-as.data.frame(cbind(hnew,llgrid0=llfit,beta1.0=betafit))
  # index2<-as.numeric(newfits$id[newfits$llgrid0==max(newfits$llgrid0)])
  
  iterpath<-suppressWarnings(dplyr::bind_rows(hpull,newfits))
  iterpath<-iterpath[order(iterpath$id),]
  rownames(iterpath)<-NULL
  if(length(printer)>0) {
    addtxt<-ifelse(length(txtadd)>0,paste("Sim:",paste(txtadd,c(", Step:"))),c("Step:"))          
    print(paste(addtxt,numiter)) 
    print(iterpath)
  }
  iterout<-rbind(iterout,data.frame(iterpath,step=numiter))
  
  index1<-as.numeric(iterpath$id[which.max(iterpath$llgrid0)])
  
  if(length(index1)>1){
    index<-ifelse(any(index1==5),5,ifelse(all(index1>=5),min(index1),max(index1)))
    rm(index1)
    index1<-index
    rm(index)
  }
  
  if(index1!=5){
    diff.llike<-abs(iterpath$llgrid0[iterpath$id==index1]-ll1)
  }else if(index1==5){
    iternewpath<-iterpath[iterpath$id!=5,]
    llike2<-max(iternewpath$llgrid0)
    diff.llike<-abs(llike2-ll1)
    rm(llike2,iternewpath)
  }
  
  diff.inhf<-abs(iterpath$h.in[iterpath$id==index1]-hin1)
  diff.outhf<-abs(iterpath$h.out[iterpath$id==index1]-hout1)
  ll1<-iterpath$llgrid0[iterpath$id==index1]
  beta1<-iterpath$beta1.0[iterpath$id==index1]
  hin1<-iterpath$h.in[iterpath$id==index1]
  hout1<-iterpath$h.out[iterpath$id==index1]
  cent<-paste(round(iterpath$h.in[iterpath$id==5],2),round(iterpath$h.out[iterpath$id==5],2),sep="|")
  
  llinm1<-which(hnew.in==hin1)
  lloutm1<-which(hnew.out==hout1)
  llinm2<-0
  lloutm2<-0
  
  while(diff.llike > tol.ll & numiter < 100 & llinm1!=llinm2 & lloutm1!=lloutm2){
    numiter<-numiter+1
    if(length(iterh)>0){
      # Output the initialized half-life values
      iter.hist<-rbind(iter.hist,data.frame(diff.inhf=round(diff.inhf,2),
                                            diff.outhf=round(diff.outhf,2),
                                            half.in=round(hin1,2),
                                            half.out=round(hout1,2),
                                            llike=round(ll1,2),
                                            diff.llike=round(diff.llike,3),
                                            beta1=round(beta1,3),
                                            Ratio=round(exp(beta1*stdose),3),
                                            center=cent))
    }
    
    # Create vector of values for fixed decline
    iter.in<-iterpath[iterpath$h.out==hout1,]
    
    # Find max and second max for incline param at fixed decline
    llinm1<-as.numeric(which.max(iter.in$llgrid0)[1])
    llinm2<-as.numeric(which(iter.in$llgrid0==max(iter.in$llgrid0[-llinm1]))[1])
    h1.in<-hin1
    h2.in<-iter.in$h.in[llinm2]
    
    # Set new incline value to the mean of the two maxed values
    h3.in<-mean(c(h1.in,h2.in))
    
    # Create new incline vector that includes two maxes bracketing mean for fixed decline
    hnew.in<-sort(c(h1.in,h3.in,h2.in))
    
    # Create vector of values for fixed incline
    iter.out<-iterpath[iterpath$h.in==hin1,]
    
    # Find max and second max for decline param at fixed incline
    lloutm1<-as.numeric(which.max(iter.out$llgrid0)[1])
    lloutm2<-as.numeric(which(iter.out$llgrid0==max(iter.out$llgrid0[-lloutm1]))[1])
    h1.out<-hout1
    h2.out<-iter.out$h.out[lloutm2]
    
    # Set new decline value to the mean of the two maxed values
    h3.out<-mean(c(h1.out,h2.out))
    
    # Create new decline vector that includes two maxes bracketing mean for fixed incline
    hnew.out<-sort(c(h1.out,h3.out,h2.out))
    
    # Setup a new grid with the refined incline and decline parameters
    hnew<-expand.grid(h.in=hnew.in,h.out=hnew.out)
    
    # Name the combinations to minimize calculations
    hnew$hname<-paste(hnew$h.in,hnew$h.out,sep="|")
    hnew$id<-rownames(hnew)
    
    iterpath<-iterpath[,names(iterpath)!=c("id")]
    
    tout<-intersect(hnew$hname,iterpath$hname)
    hpull<-merge(hnew,iterpath,by=c("hname","h.in","h.out")) # find way to merge in the logL, h.in, h.out, and beta param only
    hnew<-hnew[!(hnew$hname %in% tout),] # remove the repeated obs
    
    newdim<-length(hnew$h.in)
    
    # Fit the new models
    Cfit<-sapply(1:newdim,function(i) C1fun.2h(thalf=c(hnew$h.in[i],hnew$h.out[i]),dat=datter)$Conc)
    
    if(modtype=="cph") {
      mfit<-lapply(1:newdim,function(i) quick.cph(C.in=Cfit[,i],dat=datter,t1=t1,t2=t2,case=case,strats=s,covs=cv))
      llfit<-sapply(1:newdim,function(i) mfit[[i]]$logLik)
      betafit<-sapply(1:newdim,function(i) mfit[[i]]$beta.C)
    } else {
      mfit<-lapply(1:newdim,function(i) fit.pool(C.in=Cfit[,i],dat=datter,covs=covar))
      llfit<-sapply(1:newdim,function(i) as.numeric(logLik(mfit[[i]])))
      betafit<-sapply(1:newdim,function(i) as.numeric(mfit[[i]]$coefficients["C.in"]))
    }
    
    rm(Cfit,mfit)
    
    newfits<-as.data.frame(cbind(hnew,llgrid0=llfit,beta1.0=betafit))
    
    iterpath<-suppressWarnings(dplyr::bind_rows(hpull,newfits))
    iterpath<-iterpath[order(iterpath$id),]
    rownames(iterpath)<-NULL
    if(length(printer)>0) {
      addtxt<-ifelse(length(txtadd)>0,paste("Sim:",paste(txtadd,c(", Step:"))),c("Step:"))          
      print(paste(addtxt,numiter)) 
      print(iterpath)
    }
    iterout<-rbind(iterout,data.frame(iterpath,step=numiter))
    
    index1<-as.numeric(iterpath$id[which.max(iterpath$llgrid0)])
    
    if(length(index1)>1){
      index<-ifelse(any(index1==5),5,ifelse(all(index1>=5),min(index1),max(index1)))
      rm(index1)
      index1<-index
      rm(index)
    }
    
    adjust<-ifelse(max(iterpath$llgrid0)==ll1,1,0)
    if(adjust==0){
      diff.llike<-abs(iterpath$llgrid0[iterpath$id==index1]-ll1)
    }else if(adjust==1){
      iternewpath<-iterpath[iterpath$id!=index1,]
      llike2<-max(iternewpath$llgrid0)
      diff.llike<-abs(llike2-ll1)
      rm(llike2,iternewpath)
    }
    rm(adjust)
    
    diff.inhf<-abs(iterpath$h.in[iterpath$id==index1]-hin1)
    diff.outhf<-abs(iterpath$h.out[iterpath$id==index1]-hout1)
    ll1<-iterpath$llgrid0[iterpath$id==index1]
    beta1<-iterpath$beta1.0[iterpath$id==index1]
    hin1<-iterpath$h.in[iterpath$id==index1]
    hout1<-iterpath$h.out[iterpath$id==index1]
    cent<-paste(round(iterpath$h.in[iterpath$id==5],2),round(iterpath$h.out[iterpath$id==5],2),sep="|")
  }
  
  half.in<-iterpath$h.in[iterpath$id==index1]
  half.out<-iterpath$h.out[iterpath$id==index1]
  fin.h<-cbind(half.in,half.out)
  
  # Fit the final models
  Cfin<-C1fun.2h(thalf=c(half.in,half.out),dat=datter)$Conc
  
  if(modtype=="cph") {
    if(!is.null(stratas)) {
      stratnew<-sapply(1:length(stratas),function(t) paste0("strata(",paste0(stratas[t],")")))
      covar2<-c(covar,stratnew)
    } else {
      covar2<-covar
    }
    m.1<-fit.cph(C.in=Cfin,dat=datter,covs=covar2)
  } else {
    m.1<-fit.pool(C.in=Cfin,dat=datter,covs=covar)
  }
  
  llike<-round(as.numeric(logLik(m.1)),4)
  AIC<-2*(length(m.1$coefficients)+2)-2*llike # AIC needs additional parameters for two half-life estimates
  beta1<-as.numeric(m.1$coefficients["C.in"])
  
  beta1.se<-ifelse(modtype=="cph",
                   summary(m.1)$coefficients["C.in","se(coef)"],
                   summary(m.1)$coefficients["C.in","Std. Error"])
  
  Ratio<-round(exp(beta1*stdose),4)
  Ratio.95L<-round(exp(beta1*stdose-1.96*beta1.se*stdose),4)
  Ratio.95U<-round(exp(beta1*stdose+1.96*beta1.se*stdose),4)
  
  out <- list()
  if(length(iterh)>0){ 
    out$iter.hist<-iter.hist
    out$niter<-dim(iter.hist)[1]
  }
  
  out$final.results<-c(Beta1.Est=round(beta1,4),
                       Beta1.SE=round(beta1.se,4),
                       Ratio=Ratio,
                       Ratio.95L=Ratio.95L,
                       Ratio.95U=Ratio.95U,
                       Half.in=round(half.in,2),
                       Half.out=round(half.out,2),
                       h0=h.0,
                       logL=llike,
                       AIC=round(AIC,4),
                       numiter=numiter,
                       modtype=modtype,
                       estalgo="2param",
                       RatioDose=stdose)
  out$final.mod<-summary(m.1)
  out$extra<-iterout
  
  return(out)
}