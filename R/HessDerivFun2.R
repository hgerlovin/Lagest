#' Hessian-derivation functions
#' 
#' Functions here are for calculating the first and second derivatives for the 1 and 2 lag parameter Effective Exposure. Used in the half-life finding algorithms to compute SE. Continuation from the initial Hessian-derivation functions.
#' @param b1 Final beta parameter estimate using the estimated half-life (or half-lives)
#' @param Ct Vector or column of dataframe with the calculated effective exposure values
#' @param C2,dC.1,dC.2 Vector or column of dataframe with the calculated first derivative of the effective exposure. C2 is used in the OPEE framework, while dC.1 and dC.2 correspond to the first derivatives of the TPEE model framework effective exposure with respect to the incline and decline half-lives, respectively.
#' @param C3 Vector or column of dataframe with the calculated second derivative of the effective exposure. C3 is only applicable in the OPEE framework.
#' @param d2C.1,d2C.2,d2C.12 Vector or columns of the dataframe with the calculated second derivatives of the effective exposure values with respect to the incline, decline, and both parameters, respectively. Applicable to the TPEE framework only.
#' @param y A vector or column in the dataframe with the true event values. Should be the same length and correspond to the observations for the other input parameters. Takes a value of 1 for events and 0 for non-events.
#' @param p The predicted probability of event for a given subject-time-observation using the fully-specified model.
#' @return Output includes:
#' \describe{
#'    \item{u}{Score for the model fit under given parameters}
#'    \item{I}{Fisher's Information for the model fit under given parameters}
#'    \item{llik}{Log-likelihood for the model fit under given parameters}
#' }
#' @describeIn HessDerivFun2 Calculate the Hessian matrix components from the subject-time-specific equation values for the CPH OPEE model
#' @export
timesums1<-function(b1,Ct,C2,C3,y,p){
  A.1<-sum(p*y)
  risk<-exp(p)
  A.2b<-sum(Ct*y)
  A.2l<-b1*sum(C2*y)
  A.3<-b1*sum(C3*y)
  B.1<-sum(Ct*risk)
  B.2b<-sum(Ct*Ct*risk)
  B.2l<-sum(C2*risk*(1+b1*Ct))
  C.1<-sum(risk)
  C.2<-sum(b1*C2*risk)
  C.3<-sum(b1*risk*(C3+b1*C2*C2))
  m<-sum(y)
  
  dbeta<-A.2b-m*(B.1/C.1)
  dlam<-A.2l-m*(C.2/C.1)
  ddbeta<-m*(((B.1^2)-B.2b*C.1)/(C.1^2))
  ddlam<-A.3-m*((C.3*C.1)-(C.2^2))/(C.1^2)
  ddbl<-(A.2l/b1)-m*((B.2l*C.1)-(C.2*B.1))/(C.1^2)
  
  u<-c(dbeta,dlam)
  I<-matrix(c(ddbeta,ddbl,ddbl,ddlam),nrow=2)
  llik<-A.1-m*log(C.1)
  
  return(list(u,I,llik))
}

#' @describeIn HessDerivFun2 Calculate the Hessian matrix components from the subject-time-specific equation values for the CPH TPEE model
#' @export
timesums2<-function(b1,Ct,dC.1,dC.2,d2C.1,d2C.2,d2C.12,y,p){
  A.1<-sum(p*y)
  risk<-exp(p)
  A.2b<-sum(Ct*y)
  A.4<-b1*sum(dC.1*y)
  A.5<-b1*sum(dC.2*y)
  A.6<-b1*sum(d2C.1*y)
  A.7<-b1*sum(d2C.2*y)
  A.8<-sum(b1*d2C.12*y)
  B.1<-sum(Ct*risk)
  B.2b<-sum(Ct*Ct*risk)
  B.2.1<-sum(dC.1*risk*(1+b1*Ct))
  B.2.2<-sum(dC.2*risk*(1+b1*Ct))
  C.1<-sum(risk)
  C.4<-sum(b1*dC.1*risk)
  C.5<-sum(b1*dC.2*risk)
  C.6<-sum(b1*risk*(d2C.1+b1*dC.1*dC.1))
  C.7<-sum(b1*risk*(d2C.2+b1*dC.2*dC.2))
  C.8<-sum(b1*risk*(d2C.12+b1*dC.1*dC.2))
  m<-sum(y)
  
  dbeta<-A.2b-m*(B.1/C.1)
  dlam1<-A.4-m*(C.4/C.1)
  dlam2<-A.5-m*(C.5/C.1)
  ddbeta<-m*(((B.1*B.1)-(B.2b*C.1))/(C.1*C.1))
  ddlam1<-A.6-m*((C.6*C.1)-(C.4*C.4))/(C.1*C.1)
  ddlam2<-A.7-m*((C.7*C.1)-(C.5*C.5))/(C.1*C.1)
  ddlam12<-A.8-m*(((C.8*C.1)-(C.4*C.5))/(C.1*C.1))
  ddbl1<-(A.4/b1)-m*((B.2.1*C.1)-(C.4*B.1))/(C.1*C.1)
  ddbl2<-(A.5/b1)-m*((B.2.2*C.1)-(C.5*B.1))/(C.1*C.1)
  llik<-A.1-m*log(C.1)
  
  u<-c(dbeta,dlam1,dlam2)
  I<-matrix(c(ddbeta,ddbl1,ddbl2,ddbl1,ddlam1,ddlam12,ddbl2,ddlam12,ddlam2),nrow=3)
  return(list(u,I,llik))
}

