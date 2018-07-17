#' Basic OPEE/TPEE Calculation
#'
#' Reads in dose, start and end times, and half-life (or half-lives). Calculates concentration for inputs.
#' @param d dose/concentration value DoseX
#' @param s time since start tStartX
#' @param e time since end tEndX
#' @param h half-life parameter (for use with C1fn.h() only)
#' @param h1,h2 half-life parameters for incline and decline, respectively. For use with C1.new() only.
#' @return Outputs a single value for effective exposure concentration.
#' @describeIn C1fn.h Uses the one-parameter specification.
##' @export
C1fn.h=function(d,s,e,h){
  return(d*(exp(-log(2)*e/h)-exp(-log(2)*s/h)))
}

#' @describeIn C1fn.h Allows for two differing half-lives and calculates the effective exposure concentration based on incline and decline parameters.
#' @export
C1.new<-function(d,s,e,h1,h2){
  return(d*(1-exp(-log(2)*(s-e)/h1))*exp(-log(2)*e/h2))
}

#' OPEE/TPEE Calculation across dataset subjects and regimens
#'
#' Reads in dataset and input half-life. Uses lapply and individual effective exposure calculation method to determine the total effective exposure for each observation for each regimen given in the dataset. Sums each dose for each subject at each time point. 
#' @param thalf Assumed half-life for the OPEE calculation. Default is NULL. For the two parameter model, thalf must be input as a two-value vector with the incline half-life listed first.
#' @param dat Dataset with three columns per exposure event/regimen following the naming conventions for X total regimens: Dose1--DoseX, tStart1--tStartX, tEnd1--tEndX
#' @return Outputs a single vector of effective exposure concentrations for all observations in the dataset (should match the event vector length).
#' @describeIn C1fun.h Uses the one-parameter specification to calculated a total effective exposure for a given observation at time t. Based on sum of C1fn.h().
#' @export
C1fun.h=function(thalf=NULL,dat){
  Ntimes=length(grep("^Dose",names(dat)))
  Snames=grep("^tStart",names(dat))
  Enames=grep("^tEnd",names(dat))
  Dnames=grep("^Dose",names(dat))
  
  C1=lapply(1:Ntimes,function(q) C1fn.h(d=dat[,Dnames[q]],s=dat[,Snames[q]],e=dat[,Enames[q]],h=thalf))
  Conc=Reduce("+",lapply(C1, function(s) replace(s, is.na(s),0)))
  return(Conc)
}

#' @describeIn C1fun.h Uses the two-parameter specification to calculated a total effective exposure for a given subject at time t. Uses C1.new() instead of C1fn.h().
#' @export
C1fun.2h=function(thalf=NULL,dat){
  out<-list()
  if(length(thalf)==1) {
    out$footnote<-"Assuming the same half-life for incline and decline"
    thalf<-c(thalf,thalf)
  }
  if(length(thalf)>2) return(print("ERROR: Only one half per direction allowed"))
  
  Ntimes<-length(grep("^Dose",names(dat)))
  Snames<-grep("^tStart",names(dat))
  Enames<-grep("^tEnd",names(dat))
  Dnames<-grep("^Dose",names(dat))
  
  C1<-lapply(1:Ntimes,function(q) C1.new(d=dat[,Dnames[q]],s=dat[,Snames[q]],e=dat[,Enames[q]],h1=thalf[1],h2=thalf[2]))
  out$Conc<-Reduce("+",lapply(C1, function(s) replace(s, is.na(s),0)))
  return(out)
}