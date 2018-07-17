#' Simulate multiple-scenarios dataset
#'
#' Generate simulated datasets based on multiple pre-specified scenarios. Should only use scenarios generated using the same underlying effective exposure structure (same number of parameters and values).
#' @param combo Pre-defined scenarios to use for simulating subjects. Must be input as a list of ScenSpec/ScenSpec2 output objects.
#' @param nperscen Total number of subjects to simulate per trajectory/scenario. Default is 5,000 subjects. Can be input as a vector of values with the same number of elements as listed in "combo" - i.e. nperscen=c(5000,5000,2000) would imply 5,000 subjects generated under scenarios 1 and 2 and only 2,000 subjects with scenario 3.
#' @param seedno Random number generator seed specification for simulating event vs. probability matrix. Default is 25.
#' @inherit DatScen return
#' @seealso ScenSpec, ScenSpec2, DatSpec, DatSpec2, DatScen
#' @export
MultiScen <- function(combo=list(scen1,scen2,scen3), nperscen=5000, seedno=25) {
  scen.num<-length(combo)
  if (scen.num==1) return(print("ERROR: Only one scenario specified, use DatScen() function"))  
  
  if(length(nperscen)==1) nperscen<-rep(nperscen,scen.num)
  if(length(nperscen)!=scen.num) return(print("ERROR: NPERSCEN must have length of 1 or equal to number of scenarios"))
  
  set.seed(seedno)
  
  rs<-lapply(1:scen.num, function(s) matrix(runif(length(combo[[s]]$time)*nperscen[s]),nrow=nperscen[s],ncol=length(combo[[s]]$time)))
  probs<-lapply(1:scen.num, function(sn) combo[[sn]][,"prob"])
  
  newerdat<-list()
  events<-list()
  lt<-sapply(1:scen.num, function(s) as.numeric(length(combo[[s]]$time)))
  for(snum in 1:scen.num){
    events[[snum]]<-t(apply(sweep(rs[[snum]],2,probs[[snum]]), 1, 
                            function(s) { a=which(s<0)
                            if(length(a)<1) {
                              rep(FALSE, lt[snum])
                            } else if(min(a)==1) {
                              c(TRUE,rep(NA, lt[snum]-1))
                            } else if(min(a)==lt[snum]) {
                              c(rep(FALSE, lt[snum]-1),TRUE)
                            } else {
                              c(rep(FALSE, min(a)-1),TRUE,rep(NA, lt[snum]-min(a)))
                            }
                            }))
    
    newerdat[[snum]]<-reshape::melt.matrix(events[[snum]])
    
    colnames(newerdat[[snum]])<-c("ID","time","event")
    newerdat[[snum]]<-newerdat[[snum]][!is.na(newerdat[[snum]]$event),]
    newerdat[[snum]]<-merge(newerdat[[snum]],combo[[snum]][,!(names(combo[[snum]]) %in% c("half","prob","OR","baser"))],by=c("time"))
    newerdat[[snum]]<-newerdat[[snum]][with(newerdat[[snum]], order(ID,time)),]
    rownames(newerdat[[snum]])<-NULL
    
    newerdat[[snum]]$scen.num<-snum
  }
  
  rm(rs,events,probs,lt)
  IDvec<-sapply(1:scen.num, function(f) max(newerdat[[f]]$ID))
  for(snum in 2:scen.num){
    newerdat[[snum]]$ID<-newerdat[[snum]]$ID+sum(IDvec[1:(snum-1)])
  }
  
  fin.dat<-do.call(plyr::rbind.fill,newerdat)
  rownames(fin.dat)=NULL
  
  return(fin.dat)
}

#' @describeIn MultiScen Used for parallel processing in simulations to make a survival-based dataset following the MultiScen() function. This step removes all time-point observations that do not have a corresponding event - i.e. Only unique event-times are retained. This combines the MultiScen() and full.surv() functions to output a single dataset for simulation purposes.
#' @export
makeDat.cph<-function(combo=list(base1.up,base1.down,base1.ctrl),nperscen=5000,seedno=seedset){
  fullbase1<-MultiScen(combo=combo,nperscen=nperscen,seedno=seedno) 
  fullbase1<-fullbase1[order(fullbase1$time,fullbase1$ID),]
  ts<-unique(sort(fullbase1$time[fullbase1$event==TRUE]))
  
  studyt<-max(fullbase1$time)
  if(ts[1]==0) ts<-ts[-1]
  if(ts[length(ts)]!=studyt) ts<-c(ts,studyt)
  
  survbase1<-fullbase1[fullbase1$time %in% ts,]
  rm(fullbase1)
  
  survbase1$tstop<-survbase1$time
  survbase1<-survbase1[,!(names(survbase1) %in% c("time","intlen"))]
  survbase1<-survbase1[order(survbase1$tstop,survbase1$ID),]
  
  loc<-match(ts,survbase1$tstop)
  loc1<-loc[-1]
  loc2<-c(loc[-c(1,2)]-1,length(survbase1$tstop))
  rm(loc)
  repl<-loc2-loc1+1
  tnew<-rep(ts[-length(ts)], times=repl)
  repz<-as.numeric(length(survbase1$tstop)-length(tnew))
  survbase1$time<-c(rep(0,repz),tnew)
  rm(tnew,repz,ts,loc1,loc2,repl,studyt)
  
  survbase1<-survbase1[!is.na(survbase1$tstop),]
  
  survbase1<-survbase1[order(survbase1$time,survbase1$ID),]
  return(survbase1)
}