#' Fit Cox Proportional Hazards Regression or Pooled Logistic Regression Models
#' 
#' Basic models fit using survival and speedglm packages. For the survival models, Breslow handling of ties is used.
#' @param C.in Vector of values for the effective exposure that corresponds in length and order with the either Surv(time,tstop,event) outcome matrix (fit.cph) or length of the event column in the dataframe \code{dat}.
#' @param dat Dataframe with appropriately named columns to be used in either "Surv(time,tstop,event)~" or "event~" and that contains all covariate columns listed in the \code{covs} input. For the fit.cph(), it is assumed that the dataset is already in the appropriate survival dataset format. 
#' @param covs Vector of names of the covariates to use in the model. Default is none. Function collapses the vector into a formula object separated by "+". For fit.cph() - To specify stratifying or cluster parameters, just include the names in the vector as strata() or cluster() - i.e. covs=c("bmi","strata(age)","hypertension"). Within the find.1half() and find.2half() algorithms, the appropriate strata should be converted to work with this input.
#' @return Outputs a model object corresponding to the function used. 
#' @rdname ModFits
#' @export
fit.cph<-function(C.in,dat,covs=NULL){
  if(is.null(covs)) {
    mod<-survival::coxph(as.formula("Surv(time,tstop,event) ~ C.in"),data=data.frame(cbind(dat,C.in)),ties="breslow")
  } else {
    mod<-survival::coxph(as.formula(paste("Surv(time,tstop,event) ~ C.in + ",paste(covs,collapse="+"))),data=data.frame(cbind(dat,C.in)),ties="breslow")
  }
  return(mod)
}

#' @rdname ModFits
#' @export
fit.pool<-function(C.in,dat,covs=NULL){
  if(is.null(covs)) {
    mod<-speedglm::speedglm(as.formula("event ~ C.in"),family=binomial(),data=data.frame(cbind(dat,C.in)))
  } else {
    mod<-speedglm::speedglm(as.formula(paste("event ~ C.in + ",paste(covs,collapse="+"))),family=binomial(),data=data.frame(cbind(dat,C.in)))
  }
  return(mod)
}
