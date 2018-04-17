#' Basic OPEE Calculation
#'
#' Reads in dose, start and end times, and half-life. Calculates concentration for inputs.
#' @param d dose/concentration value DoseX
#' @param s time since start tStartX
#' @param e time since end tEndX
#' @param h half-life parameter
#' @return Outputs a numeric vector of values for concentration that is equal in length to the number of rows/observations in the dataframe \code{dat}.
#' @export
#' @examples
#' C1fn.h()

C1fn.h=function(d,s,e,h){
  return(d*(exp(-log(2)*e/h)-exp(-log(2)*s/h)))
}

#' Basic TPEE Calculation
#'
#' Reads in dose, start and end times, and two half-life values - one in, one out. Calculates concentration for inputs. Allows for two differing half-lives.
#' @param d dose/concentration value DoseX
#' @param s time since start tStartX
#' @param e time since end tEndX
#' @param h1 half-life parameter for incline
#' @param h2 half-life parameter for decline
#' @export
#' @examples
#' C1.new()

C1.new<-function(d,s,e,h1,h2){
  return(d*(1-exp(-log(2)*(s-e)/h1))*exp(-log(2)*e/h2))
}

#' OPEE Calculation across dataset subjects and regimens
#'
#' Reads in dataset and input half-life. Uses lapply and C1fn.h to calculate concentration for each observation for each regimen given in the dataset. Sums each dose for each subject at each time point. Outputs vector of concentrations for all observations in the dataset (should match the event vector length).
#' @param thalf Assumed half-life for the OPEE calculation. Default is NULL.
#' @param dat Dataset with three columns per exposure event/regimen following the naming conventions for X total regimens: Dose1--DoseX, tStart1--tStartX, tEnd1--tEndX
#' @export
#' @examples
#' C1fun.h()

C1fun.h=function(thalf=NULL,dat){
  Ntimes=length(grep("^Dose",names(dat)))
  Snames=grep("^tStart",names(dat))
  Enames=grep("^tEnd",names(dat))
  Dnames=grep("^Dose",names(dat))
  
  C1=lapply(1:Ntimes,function(q) C1fn.h(d=dat[,Dnames[q]],s=dat[,Snames[q]],e=dat[,Enames[q]],h=thalf))
  Conc=Reduce("+",lapply(C1, function(s) replace(s, is.na(s),0)))
  return(Conc)
}

#' TPEE Calculation across dataset subjects and regimens
#'
#' Reads in dataset and input half-life. Uses lapply and C1fn.h to calculate concentration for each observation for each regimen given in the dataset. Sums each dose for each subject at each time point. Outputs vector of concentrations for all observations in the dataset (should match the event vector length).
#' @param thalf Assumed pair of half-lives for the TPEE calculation. Default is NULL. Entering a single half-life results in a warning and OPEE calculation.
#' @param dat Dataset with three columns per exposure event/regimen following the naming conventions for X total regimens: Dose1--DoseX, tStart1--tStartX, tEnd1--tEndX
#' @export
#' @examples
#' C1fun.2h()

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