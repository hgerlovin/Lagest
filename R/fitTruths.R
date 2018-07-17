#' Parallel processing fits for true exposure parameters
#' 
#' Simulation-based model fitting for "true" exposure parameters.
#' @param datin Dataframe that consists of columns with the following specifications. Should be organized as one observation per subject per time, and should be in a survival dataset format for using fitTruth.cph().
#' @param a Pointer to the dataframe column with the true effective exposure for each individual at each time.
#' @param b Pointer to the dataframe column that provides the current dosing level. This is similar to the "currently" exposed model, without adjustment for past exposure, when the dose is binary.
#' @param c Pointer to the dataframe column that indicates "Ever" exposed
#' @param up Pointer to the dataframe column that indicates "Currently" exposed
#' @param dn Pointer to the dataframe column that indicates "Past" exposed
#' @param t1 Pointer to the dataframe column with "time".
#' @param t2 Pointer to the dataframe column with "tstop". Only relevant for the fitTruth.cph() function.
#' @param case Pointer to the dataframe column with "event"
#' @param truep Specifies the number of "lag" parameters to assign when calculating the AIC for the true models. Default is 1. *Look into this for the cox model fits*
#' 
#' @return Results in a dataframe consisting of model results for the true parameter fits.
#' \describe{
#'    \item{Beta1.Est and Beta1.SE}{Estimated effect size or coefficient on the log-scale, with corresponding variability.}
#'    \item{Ratio, Ratio.95L, and Ratio.95U}{Resulting estimates on the Odds or Hazard Ratio scale.}
#'    \item{logL and AIC}{The corresponding model fit values}
#'    \item{modtype}{Whether pooled logistic regression (pool) or cox proportional hazards (cph) was used}
#'    \item{estalgo}{Specification for the true parameter fit}
#'    \item{Half.in, Half.out, h0, numiter, and RatioDose}{Additional columns provided in the output to be able to merge with algorithmic results. All should be "NA"}
#' }
#' 
#' @describeIn fitTruths Fits CoxPH models for the "truth" parameters.
#' @export
fitTruth.cph<-function(datin,a,b,c,up,dn,t1,t2,case){
  
  time.s<-proc.time()
  true.res<-survival::agreg.fit(x=as.matrix(datin[,a]), y=Surv(datin[,t1],datin[,t2],datin[,case]), strata=NULL, method=1, offset=NULL, init=NULL, weights=NULL, control=coxph.control(), rownames=NULL)
  time.t<-round((proc.time()-time.s)[3],1)
  fittrue<-list()
  fittrue$logLik<-true.res$loglik[2]
  fittrue$beta.C<-as.numeric(true.res$coefficients[1])
  fittrue$beta.C.se<-sqrt(true.res$var[1,1])
  attr(fittrue,"class")[1]="cph"
  rm(true.res)
  new.true<-tryCatch({as.data.frame(cbind(getSummary(fittrue,dose=1),estalgo="trueEE",fit.time=time.t))},  error=function(e) as.data.frame(cbind(Beta1.Est=NA,Beta1.SE=NA,Ratio=NA,Ratio.95L=NA,Ratio.95U=NA,Half.in=NA,Half.out=NA,h0=NA,logL=NA,AIC=NA,numiter=NA,modtype="cph",estalgo="trueEE",RatioDose=1)))
  rm(fittrue,time.s,time.t)
  
  time.s<-proc.time()
  true.dose<-survival::agreg.fit(x=as.matrix(datin[,b]), y=Surv(datin[,t1],datin[,t2],datin[,case]), strata=NULL, method=1, offset=NULL, init=NULL, weights=NULL, control=coxph.control(), rownames=NULL)
  time.t<-round((proc.time()-time.s)[3],1)
  fitdose<-list()
  fitdose$logLik<-true.dose$loglik[2]
  fitdose$beta.C<-as.numeric(true.dose$coefficients[1])
  fitdose$beta.C.se<-sqrt(true.dose$var[1,1])
  attr(fitdose,"class")[1]="cph"
  rm(true.dose)
  new.dose<-tryCatch({as.data.frame(cbind(getSummary(fitdose,dose=1),estalgo="DoseMod",fit.time=time.t))}, error=function(e) as.data.frame(cbind(Beta1.Est=NA,Beta1.SE=NA,Ratio=NA,Ratio.95L=NA,Ratio.95U=NA,Half.in=NA,Half.out=NA,h0=NA,logL=NA,AIC=NA,numiter=NA,modtype="cph",estalgo="DoseMod",RatioDose=1)))
  rm(fitdose,time.s,time.t)
  
  time.s<-proc.time()
  true.ever<-survival::agreg.fit(x=as.matrix(datin[,c]), y=Surv(datin[,t1],datin[,t2],datin[,case]), strata=NULL, method=1, offset=NULL, init=NULL, weights=NULL, control=coxph.control(), rownames=NULL)
  time.t<-round((proc.time()-time.s)[3],1)
  fitever<-list()
  fitever$logLik<-true.ever$loglik[2]
  fitever$beta.C<-as.numeric(true.ever$coefficients[1])
  fitever$beta.C.se<-sqrt(true.ever$var[1,1])
  attr(fitever,"class")[1]="cph"
  rm(true.ever)
  new.ever<-tryCatch({as.data.frame(cbind(getSummary(fitever,dose=1),estalgo="EverMod",fit.time=time.t))}, error=function(e) as.data.frame(cbind(Beta1.Est=NA,Beta1.SE=NA,Ratio=NA,Ratio.95L=NA,Ratio.95U=NA,Half.in=NA,Half.out=NA,h0=NA,logL=NA,AIC=NA,numiter=NA,modtype="cph",estalgo="EverMod",RatioDose=1)))
  rm(fitever,time.s,time.t)
  
  time.s<-proc.time()
  true.cpn<-survival::agreg.fit(x=as.matrix(datin[,c(up,dn)]), y=Surv(datin[,t1],datin[,t2],datin[,case]), strata=NULL, method=1, offset=NULL, init=NULL, weights=NULL, control=coxph.control(), rownames=NULL)
  time.t<-round((proc.time()-time.s)[3],1)
  fitcpn<-list()
  fitcpn$logLik<-true.cpn$loglik[2]
  fitcpn2<-list()
  fitcpn2$logLik<-true.cpn$loglik[2]
  fitcpn$beta.C<-as.numeric(true.cpn$coefficients[1])
  fitcpn$beta.C.se<-sqrt(true.cpn$var[1,1])
  fitcpn2$beta.C<-as.numeric(true.cpn$coefficients[2])
  fitcpn2$beta.C.se<-sqrt(true.cpn$var[2,2])
  attr(fitcpn,"class")[1]="cph"
  attr(fitcpn2,"class")[1]="cph"
  rm(true.cpn)
  new.cpn1<-tryCatch({as.data.frame(cbind(getSummary(fitcpn,dose=1),estalgo="CatModUp",fit.time=time.t))}, error=function(e) as.data.frame(cbind(Beta1.Est=NA,Beta1.SE=NA,Ratio=NA,Ratio.95L=NA,Ratio.95U=NA,Half.in=NA,Half.out=NA,h0=NA,logL=NA,AIC=NA,numiter=NA,modtype="cph",estalgo="CatModUp",RatioDose=1)))
  new.cpn2<-tryCatch({as.data.frame(cbind(getSummary(fitcpn2,dose=1),estalgo="CatModDown",fit.time=time.t))}, error=function(e) as.data.frame(cbind(Beta1.Est=NA,Beta1.SE=NA,Ratio=NA,Ratio.95L=NA,Ratio.95U=NA,Half.in=NA,Half.out=NA,h0=NA,logL=NA,AIC=NA,numiter=NA,modtype="cph",estalgo="CatModDown",RatioDose=1)))
  rm(fitcpn,fitcpn2,time.s,time.t)
  
  new.cpn<-rbind(new.cpn1,new.cpn2)  
  fin<-suppressWarnings(dplyr::bind_rows(list(new.true,new.dose,new.ever,new.cpn)))
  fin
}

#' @describeIn fitTruths Fits Pooled logistic regression models for the "truth" parameters.
#' @export
fitTruth.pl<-function(datin,a,b,c,up,dn,case,truep=1){
  
  time.s<-proc.time()
  true.res<-speedglm::speedglm(datin[,c(case)]~datin[,c(a)],family=binomial(),data=datin)
  time.t<-round((proc.time()-time.s)[3],1)
  new.true<-tryCatch({as.data.frame(cbind(getSummary(true.res,dose=1),estalgo="trueEE",fit.time=time.t))},  error=function(e) as.data.frame(cbind(Beta1.Est=NA,Beta1.SE=NA,Ratio=NA,Ratio.95L=NA,Ratio.95U=NA,logL=NA,modtype="pool",estalgo="trueEE",RatioDose=1)))
  llike<-as.numeric(as.character(new.true$logL))
  new.true$AIC<-2*(2+truep)-2*llike # Intercept, Beta and number of lag param in true model
  rm(time.s,time.t,llike)
  
  time.s<-proc.time()
  true.dose<-speedglm::speedglm(datin[,c(case)]~datin[,c(b)],family=binomial(),data=datin)
  time.t<-round((proc.time()-time.s)[3],1)
  new.dose<-tryCatch({as.data.frame(cbind(getSummary(true.dose,dose=1),estalgo="DoseMod",fit.time=time.t))}, error=function(e) as.data.frame(cbind(Beta1.Est=NA,Beta1.SE=NA,Ratio=NA,Ratio.95L=NA,Ratio.95U=NA,logL=NA,modtype="pool",estalgo="DoseMod",RatioDose=1)))
  llike<-as.numeric(as.character(new.dose$logL))
  new.dose$AIC<-4-2*llike # Intercept and Beta estimated
  rm(llike,time.s,time.t)
  
  time.s<-proc.time()
  true.ever<-speedglm::speedglm(datin[,c(case)]~datin[,c(c)],family=binomial(),data=datin)
  time.t<-round((proc.time()-time.s)[3],1)
  new.ever<-tryCatch({as.data.frame(cbind(getSummary(true.ever,dose=1),estalgo="EverMod",fit.time=time.t))}, error=function(e) as.data.frame(cbind(Beta1.Est=NA,Beta1.SE=NA,Ratio=NA,Ratio.95L=NA,Ratio.95U=NA,logL=NA,modtype="pool",estalgo="EverMod",RatioDose=1)))
  llike<-as.numeric(as.character(new.ever$logL))
  new.ever$AIC<-4-2*llike # Intercept and Beta estimated
  rm(time.s,time.t,llike)
  
  time.s<-proc.time()
  true.cpn<-speedglm::speedglm(datin[,c(case)]~datin[,c(up)]+datin[,c(dn)],family=binomial(),data=datin)
  time.t<-round((proc.time()-time.s)[3],1)
  
  llike=as.numeric(logLik(true.cpn)[1])
  AIC=6-2*llike # Two estimated parameters in the model + intercept
  b1.c=as.numeric(summary(true.cpn)$coefficients[2,"Estimate"])
  OR.c=as.numeric(exp(b1.c))
  
  b1.c.se=as.numeric(summary(true.cpn)$coefficients[2,"Std. Error"])
  
  OR.c.95L=as.numeric(exp(b1.c-1.96*b1.c.se))
  OR.c.95U=as.numeric(exp(b1.c+1.96*b1.c.se))
  
  out.c<-cbind(round(b1.c,4),round(b1.c.se,4),round(OR.c,4),round(OR.c.95L,4),
               round(OR.c.95U,4),round(llike,4),round(AIC,4))
  
  b1.p=as.numeric(summary(true.cpn)$coefficients[3,"Estimate"])
  OR.p=as.numeric(exp(b1.p))
  
  b1.p.se=as.numeric(summary(true.cpn)$coefficients[3,"Std. Error"])
  
  OR.p.95L=as.numeric(exp(b1.p-1.96*b1.p.se))
  OR.p.95U=as.numeric(exp(b1.p+1.96*b1.p.se))
  
  out.p<-cbind(round(b1.p,4),round(b1.p.se,4),round(OR.p,4),round(OR.p.95L,4),
               round(OR.p.95U,4),round(llike,4),round(AIC,4))
  
  outnames=c("Beta1.Est","Beta1.SE","Ratio","Ratio.95L","Ratio.95U","logL","AIC")
  colnames(out.c)<-outnames
  colnames(out.p)<-outnames
  
  new.cpn1<-as.data.frame(cbind(out.c,modtype="pool",estalgo="CatModUp",fit.time=time.t,RatioDose=1))
  new.cpn2<-as.data.frame(cbind(out.p,modtype="pool",estalgo="CatModDown",fit.time=time.t,RatioDose=1))
  rm(time.s,time.t)
  
  new.cpn<-rbind(new.cpn1,new.cpn2)  
  new.cpn$AIC<-as.numeric(as.character(new.cpn$AIC))
  fin<-suppressWarnings(dplyr::bind_rows(new.true,new.dose,new.ever,new.cpn))
  fin
}

#' @describeIn fitTruths Pulls the relevant columns for analyzing "truth" values from a survival dataset
#' @export
pullTruth.cph<-function(datin){
  datout<-datin[,names(datin) %in% c("time","tstop","event","currC","currD","everD","scen.num","ID")]
  return(datout)
}

