#' Standard error approximations for Cox Proportional Hazard Model OPEE/TPEE results
#'
#' Calculates the Hessian-Derived standard errors for Beta, and Half-Life parameters (single half-life or incline/decline half-lives). Outputs a vector of values (see output descriptions). 
#' @param datter dataset with columns for the survival outcomes (time, tstop, event), effective exposure components (Dose1 -- DoseX, tStart1 -- tStartX, tEnd1 -- tEndX), covariates and stratas
#' @param t1,t2  column pointers for "time" and "tstop" in the \code{datter} dataframe
#' @param case column pointer for "event" in the \code{datter} dataframe
#' @param cv column pointers for all model covariates in the datter dataframe
#' @param s column pointers for all model strata in the datter dataframe. Only applicable for survival models
#' @param h,hin,hout Final derived estimated value for half-life (h) or incline/decline half-lives (hin, hout) parameters for which standard errors need to be approximated.
#' @param stdose Relative steady-state dose for the beta parameter interpretation
#' 
#' @return The following values are computed and returned
#' \describe{
#'    \item{Beta1.HSE}{Hessian-derived standard error for the beta estimate. Tends to be slightly larger than the estimated standard error from Cox regression.}
#'    \item{Ratio.H95L and Ratio.H95U}{Lower and Upper bounds, respectively, for Hazard Ratio based on 95\% normal approximation using hessian-derived standard errors.}
#'    \item{In.HSE}{Hessian-derived standard error for incline half-life parameter estimate. When using the OPEE model inputs, this represents the single parameter half-life standard error.}
#'    \item{In.95L and In.95U}{Lower and Upper bounds, respectively, for incline half-life parameter based on 95\% normal approximation using hessian-derived standard errors. When using the OPEE model inputs, this represents the single parameter half-life 95\% normally-approximated bounds.}
#'    \item{Out.HSE}{Hessian-derived standard error for decline half-life parameter estimate. Returns NA for OPEE model estimation.}
#'    \item{Out.95L and Out.95U}{Lower and Upper bounds, respectively, for decline half-life parameter based on 95\% normal approximation using hessian-derived standard errors. Returns NA for OPEE model estimation.}
#' }
#'    
#' @describeIn BetaH.SEcph For calculating the Hessian-derived standard errors following the OPEE algorithm.
#' @export
BetaH.SEcph<-function(datter,t1,t2,case,cv,s,h,stdose){
  Clarge<-CTfun.h(thalf=h,dat=datter)
  
  if(!is.null(s)) {
    fitme<-survival::agreg.fit(x=as.matrix(cbind(Clarge[,1],datter[,c(cv)])), y=Surv(datter[,t1],datter[,t2],datter[,case]), strata=as.numeric(strata(datter[,c(s)],shortlabel=TRUE)), method=1, offset=NULL, init=NULL, weights=NULL, control=coxph.control(), rownames=NULL)
  } else {
    fitme<-survival::agreg.fit(x=as.matrix(cbind(Clarge[,1],datter[,c(cv)])), y=Surv(datter[,t1],datter[,t2],datter[,case]), strata=NULL, method=1, offset=NULL, init=NULL, weights=NULL, control=coxph.control(), rownames=NULL)
  }
  
  beta1<-as.numeric(fitme$coefficients[1])
  
  pred<-fitme$linear.predictors
  
  gett<-cbind(Clarge,pred=pred,case=datter[,case],tstop=datter[,t2],strat=datter[,c(s)])
  rm(Clarge,pred,fitme)
  
  yvec<-which(colnames(gett)=="case")
  pvec<-which(colnames(gett)=="pred")
  C<-which(colnames(gett)=="Conc")
  C2<-which(colnames(gett)=="C2")
  C3<-which(colnames(gett)=="C3")
  
  if(length(s)>0) {
    s2<-which(colnames(gett)=="strat")
    gett2<-split.data.frame(gett,f=list(gett[,"tstop"],gett[,c(s2)]))
  } else {
    gett2<-split.data.frame(gett,f=list(gett[,"tstop"]))
  }
  
  rm(gett)
  flo<-lapply(1:length(gett2), function(t) unlist(timesums1(b1=beta1, Ct=gett2[[t]][,c(C)], C2=gett2[[t]][,c(C2)], C3=gett2[[t]][,c(C3)],y=gett2[[t]][,c(yvec)],p=gett2[[t]][,c(pvec)]),use.names = F))
  
  fly<-do.call(rbind,flo)
  rm(flo)
  I1<-sum(fly[,3],na.rm=T)
  I2<-sum(fly[,4],na.rm=T)
  I3<-sum(fly[,5],na.rm=T)
  I4<-sum(fly[,6],na.rm=T)
  Imat<-matrix(c(I1,I2,I3,I4),nrow=2)
  
  sehb<-sqrt(abs(diag(solve(-Imat))))
  rm(fly,Imat)
  
  In.HSE<-sehb[2]
  Beta1.HSE<-sehb[1]
  In.95U<-h+1.96*sehb[2]
  In.95L<-h-1.96*sehb[2]
  Ratio.H95U<-exp(beta1*stdose+1.96*sehb[1]*stdose)
  Ratio.H95L<-exp(beta1*stdose-1.96*sehb[1]*stdose)
  
  return(c(Beta1.HSE=round(Beta1.HSE,4),Ratio.H95L=round(Ratio.H95L,4),Ratio.H95U=round(Ratio.H95U,4),In.HSE=round(In.HSE,3),In.95L=round(In.95L,2),In.95U=round(In.95U,2)))
}

#' @describeIn BetaH.SEcph For calculating the Hessian-derived standard errors following the TPEE algorithm.
#' @export
Beta2H.SEcph<-function(datter,t1,t2,case,cv,s,hin,hout,stdose){
  Clarge<-CTfun2.h(thalf=c(hin,hout),dat=datter)
  
  if(length(s)>0) {
    fitme<-survival::agreg.fit(x=as.matrix(cbind(Clarge[,1],datter[,c(cv)])), y=Surv(datter[,t1],datter[,t2],datter[,case]), strata=as.numeric(strata(datter[,c(s)],shortlabel=TRUE)), method=1, offset=NULL, init=NULL, weights=NULL, control=coxph.control(), rownames=NULL)
  } else {
    fitme<-survival::agreg.fit(x=as.matrix(cbind(Clarge[,1],datter[,c(cv)])), y=Surv(datter[,t1],datter[,t2],datter[,case]), strata=NULL, method=1, offset=NULL, init=NULL, weights=NULL, control=coxph.control(), rownames=NULL)
  }
  
  beta1<-as.numeric(fitme$coefficients[1])
  
  pred<-fitme$linear.predictors
  
  gett<-cbind(Clarge,pred=pred,case=datter[,case],tstop=datter[,t2],strat=datter[,c(s)])
  rm(Clarge,pred,fitme)
  
  C<-which(colnames(gett)=="C")
  C.h1<-which(colnames(gett)=="C.h1")
  C.h2<-which(colnames(gett)=="C.h2")
  C2.h1<-which(colnames(gett)=="C2.h1")
  C2.h2<-which(colnames(gett)=="C2.h2")
  C2.h1h2<-which(colnames(gett)=="C2.h1h2")
  yvec<-which(colnames(gett)=="case")
  pvec<-which(colnames(gett)=="pred")
  
  if(length(s)>0) {
    s2<-which(colnames(gett)=="strat")
    gett2<-split.data.frame(gett,f=list(gett[,"tstop"],gett[,c(s2)]))
  } else {
    gett2<-split.data.frame(gett,f=list(gett[,"tstop"]))
  }
  
  rm(gett)
  
  flo<-lapply(1:length(gett2), function(t) unlist(timesums2(b1=beta1, Ct=gett2[[t]][,c(C)], dC.1=gett2[[t]][,c(C.h1)], dC.2=gett2[[t]][,c(C.h2)], d2C.1=gett2[[t]][,c(C2.h1)],d2C.2=gett2[[t]][,c(C2.h2)],d2C.12=gett2[[t]][,c(C2.h1h2)],y=gett2[[t]][,c(yvec)],p=gett2[[t]][,c(pvec)])))
  fly<-do.call(rbind,flo)
  rm(flo)
  
  I1<-sum(fly[,4],na.rm=T)
  I2<-sum(fly[,5],na.rm=T)
  I3<-sum(fly[,6],na.rm=T)
  I4<-sum(fly[,7],na.rm=T)
  I5<-sum(fly[,8],na.rm=T)
  I6<-sum(fly[,9],na.rm=T)
  I7<-sum(fly[,10],na.rm=T)
  I8<-sum(fly[,11],na.rm=T)
  I9<-sum(fly[,12],na.rm=T)
  Imat<-matrix(c(I1,I2,I3,I4,I5,I6,I7,I8,I9),nrow=3)
  
  se2hb<-sqrt(abs(diag(solve(-Imat))))
  rm(fly,Imat)
  
  Ratio.H95L<-exp(beta1*stdose-1.96*se2hb[1]*stdose)
  Ratio.H95U<-exp(beta1*stdose+1.96*se2hb[1]*stdose)
  In.95L<-hin-1.96*se2hb[2]
  In.95U<-hin+1.96*se2hb[2]
  Out.95L<-hout-1.96*se2hb[3]
  Out.95U<-hout+1.96*se2hb[3]
  Beta1.HSE<-se2hb[1]
  In.HSE<-se2hb[2]
  Out.HSE<-se2hb[3]
  
  return(c(Beta1.HSE=round(Beta1.HSE,4),Ratio.H95L=round(Ratio.H95L,4),Ratio.H95U=round(Ratio.H95U,4),In.HSE=round(In.HSE,3),In.95L=round(In.95L,2),In.95U=round(In.95U,2),Out.HSE=round(Out.HSE,3),Out.95L=round(Out.95L,2),Out.95U=round(Out.95U,2)))
}
