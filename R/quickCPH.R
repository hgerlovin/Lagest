#' Fast Cox Proportional Hazard Model fitter for streamlined simulations and algorithms
#' 
#' Fast fitter for coxph models using the agreg.fit() function and column pointers. Defaults to Breslow's handling of ties.
#' @param C.in	Vector of values for the effective exposure that corresponds in length and order with the Surv(time,tstop,event) outcome matrix.
#' @param dat	Dataframe or matrix with pointers (t1, t2, case, covs, strats) specifying the relevant columns to use in the model.
#' @param t1	Pointer to the "time" column.
#' @param t2	Pointer to the "tstop" column.
#' @param case	Pointer to the "event" column.
#' @param covs	Vector of pointer values for the covariate columns. Default is none.
#' @param strats	Vector of pointer values for the stratifying columns. Default is none.
#' @param offs	Input specification of "offset" for agreg.fit(). Default is NULL.
#' @param inits	Input specification of "init" for agreg.fit(). Default is NULL.
#' @param wts	Input specification of "weights" for agreg.fit(). Default is NULL.
#' @param rnames	Input specification of "rownames" for agreg.fit(). Default is NULL.
#' @param cntrl	Default control setup for coxph models (coxph.control()).
#' @param type	Default (1) specifies using Breslow's handling of tied events.
#' @return Outputs a model of class="cph" with the log-likelihood (element 1), the beta estimate for C.in (element 2) and the corresponding std. error (element 3). Elements in the object are not labeled for speed purposes.
#' \describe{
#'    \item{output[1]=logLik}{Model fit log-likelihood value.}
#'    \item{output[2]=beta.C}{Model fit beta estimate for effective exposure.}
#'    \item{output[3]=beta.C.se}{Model fit standard error (sqrt(var[1,1])) for the beta estimate.}
#' } 
#' @export
quick.cph<-function(C.in,dat,t1,t2,
                    case,strats=NULL,covs=NULL, 
                    offs=NULL,inits=NULL,wts=NULL, 
                    rnames=NULL,cntrl=coxph.control(),type=1){
  
  if(length(strats)==0){
    if(length(covs)==0){
      fitme<-survival::agreg.fit(x=as.matrix(C.in), y=Surv(dat[,t1],dat[,t2],dat[,case]), strata=strats, method=type, offset=offs, init=inits, weights=wts, control=cntrl, rownames=rnames)
    } else if(length(covs)>0){
      fitme<-survival::agreg.fit(x=as.matrix(cbind(C.in,dat[,c(covs)])), y=Surv(dat[,t1],dat[,t2],dat[,case]), strata=strats, method=type, offset=offs, init=inits, weights=wts, control=cntrl, rownames=rnames)
    } 
  } else if(length(strats)>0){
    if(length(covs)==0){
      fitme<-survival::agreg.fit(x=as.matrix(C.in), y=Surv(dat[,t1],dat[,t2],dat[,case]), strata=as.numeric(strata(dat[,c(strats)],shortlabel=TRUE)), method=type, offset=offs, init=inits, weights=wts, control=cntrl, rownames=rnames)
    } else if(length(covs)>0){
      fitme<-survival::agreg.fit(x=as.matrix(cbind(C.in,dat[,c(covs)])), y=Surv(dat[,t1],dat[,t2],dat[,case]), strata=as.numeric(strata(dat[,c(strats)],shortlabel=TRUE)), method=type, offset=offs, init=inits, weights=wts, control=cntrl, rownames=rnames)
    } 
  }
  
  fityou<-list()
  
  fityou$logLik<-as.numeric(fitme$loglik[2])
  fityou$beta.C<-as.numeric(fitme$coefficients[1])
  fityou$beta.C.se<-sqrt(fitme$var[1,1])
  
  return(fityou)
}
