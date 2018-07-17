#' Standard error approximations for Pooled Logistic Regression Model OPEE/TPEE results
#'
#' Calculates the Hessian-Derived standard errors for Beta, and Half-Life parameters (single half-life or incline/decline half-lives). Outputs a vector of values (see output descriptions). *Currently does not adjust for covariates*
#' @param datset Dataset with columns for the outcome (event), effective exposure components (Dose1 -- DoseX, tStart1 -- tStartX, tEnd1 -- tEndX), and covariates. 
#' @param inh,outh Final derived estimated value(s) for single half-life (inh) or incline/decline half-lives (inh,outh) parameters for which standard errors need to be approximated.
#' @param stdose Relative steady-state dose for the beta parameter interpretation. Default is binary (1) exposure.
#' 
#' @inherit BetaH.SEcph return
#'    
#' @describeIn HessPool For calculating the Hessian-derived standard errors following the OPEE algorithm.
#' @export
BetaH.SEpl<-function(inh,datset,stdose=1){
  Clarge<-CTfun.h(thalf=inh,dat=datset)
  
  fitme<-stats::glm(datset$event~Clarge[,1],family=binomial(),data=data.frame(cbind(datset,Clarge)))
  
  beta1<-as.numeric(fitme$coefficients[2])
  
  pred<-fitme$fitted.values
  
  C2<-Clarge[,2]
  C3<-Clarge[,3]
  C<-Clarge[,1]
  y<-datset$event
  ddbl<-C2*(y-pred-beta1*C*pred*(1-pred))
  ddlam<-C3*beta1*(y-pred)-(C2^2)*(beta1^2)*pred*(1-pred)
  ddbeta<-(-1)*(C^2)*pred*(1-pred)
  
  llike1<-y*log(pred)+(1-y)*log(1-pred)
  totllike<-sum(llike1,na.rm=T)
  
  I1<-sum(ddbeta,na.rm=T)
  I2<-sum(ddbl,na.rm=T)
  I3<-sum(ddbl,na.rm=T)
  I4<-sum(ddlam,na.rm=T)
  Imat<-matrix(c(I1,I2,I3,I4),nrow=2)
  sehb<-sqrt(abs(diag(solve(-Imat))))
  
  Ratio.H95L<-exp(beta1*stdose-1.96*sehb[1]*stdose)
  Ratio.H95U<-exp(beta1*stdose+1.96*sehb[1]*stdose)
  In.95L<-inh-1.96*sehb[2]
  In.95U<-inh+1.96*sehb[2]
  Beta1.HSE<-sehb[1]
  In.HSE<-sehb[2]
  
  return(c(Beta1.HSE=round(Beta1.HSE,4),Ratio.H95L=round(Ratio.H95L,4),Ratio.H95U=round(Ratio.H95U,4),In.HSE=round(In.HSE,3),In.95L=round(In.95L,2),In.95U=round(In.95U,2)))
}

#' @describeIn HessPool For calculating the Hessian-derived standard errors following the TPEE algorithm.
#' @export
Beta2H.SEpl<-function(inh,outh,datset,stdose=1){
  Clarge<-CTfun2.h(thalf=c(inh,outh),dat=datset)
  
  fitme<-stats::glm(datset$event~Clarge[,1],family=binomial(),data=data.frame(cbind(datset,Clarge)))
  
  beta1<-as.numeric(fitme$coefficients[2])
  
  pred<-fitme$fitted.values
  
  C.h1<-Clarge[,2]
  C.h2<-Clarge[,3]
  C2.h1<-Clarge[,4]
  C2.h2<-Clarge[,5]
  C2.h1h2<-Clarge[,6]
  C<-Clarge[,1]
  y<-datset$event
  
  ddbl1<-C.h1*(y-pred-beta1*C*pred*(1-pred))
  ddbl2<-C.h2*(y-pred-beta1*C*pred*(1-pred))
  ddlam1<-C2.h1*beta1*(y-pred)-(C.h1^2)*(beta1^2)*pred*(1-pred)
  ddlam2<-C2.h2*beta1*(y-pred)-(C.h2^2)*(beta1^2)*pred*(1-pred)
  ddlam12<-C2.h1h2*beta1*(y-pred)-C.h2*C.h1*(beta1^2)*pred*(1-pred)
  ddbeta<-(-1)*(C^2)*pred*(1-pred)
  
  llike1<-y*log(pred)+(1-y)*log(1-pred)
  totllike<-sum(llike1,na.rm=T)
  
  I1<-sum(ddbeta,na.rm=T)
  I2<-sum(ddbl1,na.rm=T)
  I3<-sum(ddbl2,na.rm=T)
  I4<-sum(ddbl1,na.rm=T)
  I5<-sum(ddlam1,na.rm=T)
  I6<-sum(ddlam12,na.rm=T)
  I7<-sum(ddbl2,na.rm=T)
  I8<-sum(ddlam12,na.rm=T)
  I9<-sum(ddlam2,na.rm=T)
  Imat<-matrix(c(I1,I2,I3,I4,I5,I6,I7,I8,I9),nrow=3)
  sehb<-sqrt(abs(diag(solve(-Imat))))
  
  Ratio.H95L<-exp(beta1*stdose-1.96*sehb[1]*stdose)
  Ratio.H95U<-exp(beta1*stdose+1.96*sehb[1]*stdose)
  In.95L<-inh-1.96*sehb[2]
  In.95U<-inh+1.96*sehb[2]
  Out.95L<-outh-1.96*sehb[3]
  Out.95U<-outh+1.96*sehb[3]
  Beta1.HSE<-sehb[1]
  In.HSE<-sehb[2]
  Out.HSE<-sehb[3]
  
  return(c(Beta1.HSE=round(Beta1.HSE,4),Ratio.H95L=round(Ratio.H95L,4),Ratio.H95U=round(Ratio.H95U,4),In.HSE=round(In.HSE,3),In.95L=round(In.95L,2),In.95U=round(In.95U,2),Out.HSE=round(Out.HSE,3),Out.95L=round(Out.95L,2),Out.95U=round(Out.95U,2)))
}