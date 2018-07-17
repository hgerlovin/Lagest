#' Hessian-derivation functions
#' 
#' Functions here are for calculating the first and second derivatives for the 1 and 2 lag parameter Effective Exposure. Used in the half-life finding algorithms to compute SE.
#' @inheritParams C1fn.h
#' @inheritParams C1fun.h
#' @param k,k1,k2 decay-rate parameterizations for single, incline and decline half-lives, respectively
#' @return Outputs vectors for the first and second derivatives with respect to the parameters of interest.
#' @rdname HessDerivFun
#' @export
CTfun.h=function(thalf=NULL,k=NULL,dat){
  Ntimes=length(grep("^Dose",names(dat)))
  Snames=grep("^tStart",names(dat))
  Enames=grep("^tEnd",names(dat))
  Dnames=grep("^Dose",names(dat))
  
  if(is.null(k) & is.null(thalf)) return(print("ERROR: Specify either k or thalf"))
  
  C1=lapply(1:Ntimes,function(q) C1fn.h(d=dat[,Dnames[q]],s=dat[,Snames[q]],e=dat[,Enames[q]],h=thalf))
  C1tot=Reduce("+",C1)
  c2=lapply(1:Ntimes,function(q) C2fn.h(d=dat[,Dnames[q]],s=dat[,Snames[q]],e=dat[,Enames[q]],h=thalf))
  C2tot=Reduce("+",c2)
  c3=lapply(1:Ntimes,function(q) C3fn.h(d=dat[,Dnames[q]],s=dat[,Snames[q]],e=dat[,Enames[q]],h=thalf))
  C3tot=Reduce("+",c3)
  return(cbind(Conc=C1tot,C2=C2tot,C3=C3tot))
}

#' @rdname HessDerivFun
#' @export
CTfun.k=function(thalf=NULL,k=NULL,dat){
  Ntimes=length(grep("^Dose",names(dat)))
  Snames=grep("^tStart",names(dat))
  Enames=grep("^tEnd",names(dat))
  Dnames=grep("^Dose",names(dat))
  
  if(is.null(k) & is.null(thalf)) return(print("ERROR: Specify either k or thalf"))
  
  if(is.null(k) & !is.null(thalf)) k=log(2)/thalf
  
  C1=lapply(1:Ntimes,function(q) C1fn.k(d=dat[,Dnames[q]],s=dat[,Snames[q]],e=dat[,Enames[q]],k=k))
  C1tot=Reduce("+",C1)
  c2=lapply(1:Ntimes,function(q) C2fn.k(d=dat[,Dnames[q]],s=dat[,Snames[q]],e=dat[,Enames[q]],k=k))
  C2tot=Reduce("+",c2)
  c3=lapply(1:Ntimes,function(q) C3fn.k(d=dat[,Dnames[q]],s=dat[,Snames[q]],e=dat[,Enames[q]],k=k))
  C3tot=Reduce("+",c3)
  return(cbind(Conc=C1tot,C2=C2tot,C3=C3tot))
}

#' @rdname HessDerivFun
#' @export
CTfun2.k=function(thalf=NULL,k=NULL,dat){
  Ntimes=length(grep("^Dose",names(dat)))
  Snames=grep("^tStart",names(dat))
  Enames=grep("^tEnd",names(dat))
  Dnames=grep("^Dose",names(dat))
  
  if(is.null(k) & is.null(thalf)) return(print("ERROR: Specify either k or thalf"))
  
  if(is.null(k) & !is.null(thalf)) {
    k1=log(2)/thalf[1]
    k2=log(2)/thalf[2]
  }
  
  C1=lapply(1:Ntimes,function(q) C1fn.k12(d=dat[,Dnames[q]],s=dat[,Snames[q]],e=dat[,Enames[q]],k1=k1,k2=k2))
  C1tot=Reduce("+",C1)
  
  c2.1=lapply(1:Ntimes,function(q) C2fn.k1(d=dat[,Dnames[q]],s=dat[,Snames[q]],e=dat[,Enames[q]],k1=k1,k2=k2))
  C.k1=Reduce("+",c2.1)
  c2.2=lapply(1:Ntimes,function(q) C2fn.k2(d=dat[,Dnames[q]],s=dat[,Snames[q]],e=dat[,Enames[q]],k1=k1,k2=k2))
  C.k2=Reduce("+",c2.2)
  
  c3.1=lapply(1:Ntimes,function(q) C3fn.k1(d=dat[,Dnames[q]],s=dat[,Snames[q]],e=dat[,Enames[q]],k1=k1,k2=k2))
  C2.k1=Reduce("+",c3.1)
  c3.2=lapply(1:Ntimes,function(q) C3fn.k2(d=dat[,Dnames[q]],s=dat[,Snames[q]],e=dat[,Enames[q]],k1=k1,k2=k2))
  C2.k2=Reduce("+",c3.2)
  c3.12=lapply(1:Ntimes,function(q) C3fn.k1k2(d=dat[,Dnames[q]],s=dat[,Snames[q]],e=dat[,Enames[q]],k1=k1,k2=k2))
  C2.k1k2=Reduce("+",c3.12)
  
  return(cbind(C=C1tot,C.k1=C.k1,C.k2=C.k2,C2.k1=C2.k1,C2.k2=C2.k2,C2.k1k2=C2.k1k2))
}

#' @rdname HessDerivFun
#' @export
CTfun2.h=function(thalf=NULL,k=NULL,dat){
  Ntimes=length(grep("^Dose",names(dat)))
  Snames=grep("^tStart",names(dat))
  Enames=grep("^tEnd",names(dat))
  Dnames=grep("^Dose",names(dat))
  
  if(is.null(k) & is.null(thalf)) return(print("ERROR: Specify either k or thalf"))
  
  h1=thalf[1]
  h2=thalf[2]
  
  C1=lapply(1:Ntimes,function(q) C1.new(d=dat[,Dnames[q]],s=dat[,Snames[q]],e=dat[,Enames[q]],h1=h1,h2=h2))
  C1tot=Reduce("+",C1)
  
  c2.1=lapply(1:Ntimes,function(q) C2fn.h1(d=dat[,Dnames[q]],s=dat[,Snames[q]],e=dat[,Enames[q]],h1=h1,h2=h2))
  C.h1=Reduce("+",c2.1)
  c2.2=lapply(1:Ntimes,function(q) C2fn.h2(d=dat[,Dnames[q]],s=dat[,Snames[q]],e=dat[,Enames[q]],h1=h1,h2=h2))
  C.h2=Reduce("+",c2.2)
  
  c3.1=lapply(1:Ntimes,function(q) C3fn.h1(d=dat[,Dnames[q]],s=dat[,Snames[q]],e=dat[,Enames[q]],h1=h1,h2=h2))
  C2.h1=Reduce("+",c3.1)
  c3.2=lapply(1:Ntimes,function(q) C3fn.h2(d=dat[,Dnames[q]],s=dat[,Snames[q]],e=dat[,Enames[q]],h1=h1,h2=h2))
  C2.h2=Reduce("+",c3.2)
  c3.12=lapply(1:Ntimes,function(q) C3fn.h1h2(d=dat[,Dnames[q]],s=dat[,Snames[q]],e=dat[,Enames[q]],h1=h1,h2=h2))
  C2.h1h2=Reduce("+",c3.12)
  
  return(cbind(C=C1tot,C.h1=C.h1,C.h2=C.h2,C2.h1=C2.h1,C2.h2=C2.h2,C2.h1h2=C2.h1h2))
}

#' @rdname HessDerivFun
#' @export
C1fn.k=function(d,s,e,k){
  return(d*(exp(-k*e)-exp(-k*s)))
}

#' @rdname HessDerivFun
#' @export
C1fn.k12<-function(d,s,e,k1,k2){
  return(d*(1-exp(-k1*(s-e)))*exp(-k2*e))
}

#' @rdname HessDerivFun
#' @export
C1fn.h12<-function(d,s,e,h1,h2){
  return(d*(1-exp(-log(2)*(s-e)/h1))*exp(-log(2)*e/h2))
}

#' @rdname HessDerivFun
#' @export
C2fn.k=function(d,s,e,k){
  return(d*(s*exp(-k*s)-e*exp(-k*e)))
}

#' @rdname HessDerivFun
#' @export
C2fn.h=function(d,s,e,h){
  return(d*(log(2)/(h^2))*(e*exp(-log(2)*e/h)-s*exp(-log(2)*s/h)))
}

#' @rdname HessDerivFun
#' @export
C2fn.k1<-function(d,s,e,k1,k2){
  return(d*(s-e)*exp(-k1*(s-e)-k2*e))
}

#' @rdname HessDerivFun
#' @export
C2fn.h1<-function(d,s,e,h1,h2){
  return(d*(s-e)*(-log(2)/(h1*h1))*exp(-(log(2)*(s-e)/h1))*exp(-(log(2)*e/h2)))
}

#' @rdname HessDerivFun
#' @export
C2fn.k2<-function(d,s,e,k1,k2){
  return(d*(-e)*(exp(-k2*e)-exp(-k1*(s-e)-k2*e)))
}

#' @rdname HessDerivFun
#' @export
C2fn.h2<-function(d,s,e,h1,h2){
  return(d*e*log(2)/(h2*h2)*exp(-(log(2)*e/h2))*(1-exp(-(log(2)*(s-e)/h1))))
}

#' @rdname HessDerivFun
#' @export
C3fn.k=function(d,s,e,k){
  return(d*((e^2)*exp(-k*e)-(s^2)*exp(-k*s)))
}

#' @rdname HessDerivFun
#' @export
C3fn.h=function(d,s,e,h){
  # a<-((-2)*e*log(2)/(h^3))*d*exp(-e*log(2)/h)
  # b<-((e*log(2)/(h^2))^2)*d*exp(-e*log(2)/h)
  # c<-((-2)*s*log(2)/(h^3))*d*exp(-s*log(2)/h)
  # f<-((s*log(2)/(h^2))^2)*d*exp(-s*log(2)/h)
  # return((a+b)-(c+f))
  return(d*(log(2)/(h^3))*((-2*e*exp(-log(2)*e/h))+((e^2)*(log(2)/h)*exp(-log(2)*e/h))+(2*s*exp(-log(2)*s/h))-((s^2)*(log(2)/h)*exp(-log(2)*s/h))))
}

#' @rdname HessDerivFun
#' @export
C3fn.k1<-function(d,s,e,k1,k2){
  return(-d*((s-e)^2)*exp(-k1*(s-e)-k2*e))
}

#' @rdname HessDerivFun
#' @export
C3fn.h1<-function(d,s,e,h1,h2){
  return(d*(s-e)*log(2)/(h1*h1)*exp(-log(2)*(e)/h2)*exp(-log(2)*(s-e)/h1)*((2/h1)-(log(2)*(s-e))/(h1*h1)))
}

#' @rdname HessDerivFun
#' @export
C3fn.k2<-function(d,s,e,k1,k2){
  return(d*(e^2)*(exp(-k2*e)-(exp(-k1*(s-e)-k2*e))))
}

#' @rdname HessDerivFun
#' @export
C3fn.h2<-function(d,s,e,h1,h2){
  return(d*exp(-log(2)*e/h2)*(1-exp(-log(2)*(s-e)/h1))*e*log(2)/(h2*h2)*((2/h2)-(log(2)/(h2*h2)*e)))
}

#' @rdname HessDerivFun
#' @export
C3fn.k1k2<-function(d,s,e,k1,k2){
  return(-d*e*(s-e)*exp(-k1*(s-e)-k2*e))
}

#' @rdname HessDerivFun
#' @export
C3fn.h1h2<-function(d,s,e,h1,h2){
  return(-d*(s-e)*e*(log(2))*(log(2))/(h1*h1)/(h2*h2)*exp(-log(2)*e/h2)*exp(-log(2)*(s-e)/h1))
}