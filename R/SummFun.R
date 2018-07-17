#' Functions for returning summary statistics from model objects
#' 
#' @param mod Model object. Can have the \code{class} "speedglm", "coxph", or "cph".
#' @param dose Steady-state plateau for the beta parameter interpretation
#' @return Estimated coefficients (with standard errors), exponentiated values (Ratio and 95\% CI), and log-likelihood for the model.
#' @rdname SummFun
#' @export
getSummary=function(mod,dose=1){
  out<-tryCatch({
    
    if(attr(mod,"class")[1]=="speedglm"){
      llike=as.numeric(logLik(mod)[1])
      b1=as.numeric(summary(mod)$coefficients[2,"Estimate"])
      OR=as.numeric(exp(b1*dose))
      
      b1.se=as.numeric(summary(mod)$coefficients[2,"Std. Error"])
      
      OR.95L=as.numeric(exp(b1*dose-1.96*b1.se*dose))
      OR.95U=as.numeric(exp(b1*dose+1.96*b1.se*dose))
      
      out<-cbind(round(b1,4),round(b1.se,4),round(OR,4),round(OR.95L,4),
                 round(OR.95U,4),round(llike,4),dose,"pool")
      
      outnames=c("Beta1.Est","Beta1.SE","Ratio","Ratio.95L","Ratio.95U",
                 "logL","RatioDose","modtype")
      colnames(out)<-outnames
      return(out)
    }
    
    if(attr(mod,"class")[1]=="coxph"){
      llike=as.numeric(logLik(mod)[1])
      b1<-as.numeric(mod$coefficients[1])
      HR<-as.numeric(exp(b1*dose))
      
      b1.se<-as.numeric(summary(mod)$coefficients[1,"se(coef)"])
      
      HR.95L<-as.numeric(exp(b1*dose-1.96*b1.se*dose))
      HR.95U<-as.numeric(exp(b1*dose+1.96*b1.se*dose))
      
      out<-cbind(round(b1,4),round(b1.se,4),round(HR,4),round(HR.95L,4),
                 round(HR.95U,4),round(llike,4),dose,"cph")
      
      outnames=c("Beta1.Est","Beta1.SE","Ratio","Ratio.95L","Ratio.95U",
                 "logL","RatioDose","modtype")
      colnames(out)<-outnames
      return(out)
    }
    
    if(attr(mod,"class")[1]=="cph"){
      llike=as.numeric(mod$logLik)
      b1<-as.numeric(mod$beta.C)
      HR<-as.numeric(exp(b1*dose))
      
      b1.se<-as.numeric(mod$beta.C.se)
      
      HR.95L<-as.numeric(exp(b1*dose-1.96*b1.se*dose))
      HR.95U<-as.numeric(exp(b1*dose+1.96*b1.se*dose))
      
      out<-cbind(round(b1,4),round(b1.se,4),round(HR,4),round(HR.95L,4),
                 round(HR.95U,4),round(llike,4),dose,"cph")
      
      outnames=c("Beta1.Est","Beta1.SE","Ratio","Ratio.95L","Ratio.95U",
                 "logL","RatioDose","modtype")
      colnames(out)<-outnames
      return(out)
    }
  }, error=function(e) cbind(NA, NA, NA, NA, NA, NA, NA, NA) )
  
  outnames=c("Beta1.Est","Beta1.SE","Ratio","Ratio.95L","Ratio.95U",
             "logL","RatioDose","modtype")
  colnames(out)<-outnames
  
  return(out)
}

#' @param modspec Model object to use
#' @param coefname Exact name of the variables for which Odds or Hazard ratio should be calculated. May include a vector of multiple coefficients.
#' @return returnHR() and returnOR() functions additionally output the AIC for the model fit
#' @rdname SummFun
#' @export
returnOR<-function(modspec,coefname){
  ncoef<-as.numeric(length(coefname))
  OR<-numeric(length=ncoef)
  OR.95CI<-character(length=ncoef)
  Beta<-character(length=ncoef)
  param<-character(length=ncoef)
  
  for(y in 1:ncoef){
    beta.C<-as.numeric(modspec$coefficients[coefname[y]])
    se.C<-as.numeric(summary(modspec)$coefficients[coefname[y],"Std. Error"])
    OR.95L<-round(exp(beta.C-1.96*se.C),3)
    OR.95U<-round(exp(beta.C+1.96*se.C),3)
    
    OR[y]<-round(exp(beta.C),4)
    OR.95CI[y]<-paste0("(",paste0(OR.95L,paste0("-",paste0(OR.95U,")"))))
    Beta[y]<-paste(round(beta.C,4),paste0("(",paste0(round(se.C,3),")")))
    param[y]<-coefname[y]
  }
  resout<-cbind(param=param,Beta.SE=Beta,OR=OR,OR.95CI=OR.95CI)
  modname<-deparse(substitute(modspec))
  return(cbind(resout,logL=round(logLik(modspec),2),AIC=round(AIC(modspec),2),modname))
}

#' @rdname SummFun
#' @export
returnHR<-function(modspec,coefname){
  ncoef<-as.numeric(length(coefname))
  HR<-numeric(length=ncoef)
  HR.95CI<-character(length=ncoef)
  Beta<-character(length=ncoef)
  param<-character(length=ncoef)
  
  for(y in 1:ncoef){
    beta.C<-as.numeric(modspec$coefficients[coefname[y]])
    se.C<-as.numeric(summary(modspec)$coefficients[coefname[y],"se(coef)"])
    HR.95L<-format(round(as.numeric(summary(modspec)$conf.int[coefname[y],"lower .95"]),3),nsmall=3)
    HR.95U<-format(round(as.numeric(summary(modspec)$conf.int[coefname[y],"upper .95"]),3),nsmall=3)
    
    HR[y]<-format(round(as.numeric(summary(modspec)$conf.int[coefname[y],"exp(coef)"]),4),nsmall=4)
    HR.95CI[y]<-paste0("(",paste0(HR.95L,paste0("-",paste0(HR.95U,")"))))
    Beta[y]<-paste(format(round(beta.C,4),nsmall=4),paste0("(",paste0(format(round(se.C,3),nsmall=3),")")))
    param[y]<-coefname[y]
  }
  resout<-cbind(param=param,Beta.SE=Beta,HR=HR,HR.95CI=HR.95CI)
  modname<-deparse(substitute(modspec))
  return(cbind(resout,logL=round(logLik(modspec),2),AIC=round(AIC(modspec),2),modname))
}