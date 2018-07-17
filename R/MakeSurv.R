#' Prepare dataset for survival analysis
#' 
#' @param datin Input dataframe with one observation for each subject at each study time point. Must have "time", "tstop" and "event" columns with data points for all subjects at all event times.
#' @return One observation per subject per event time. New column, tstop, takes the value of the event time, while the "time" column reflects the previous event time. I.e. Assumes right-censored times such that "tstop" values correspond to "time" from the original input dataframe and the revised "time" column reflects the preceding event time. 
#' @seealso makeDat.cph
#' @export
full.surv<-function(datin){
  ts<-unique(sort(datin$time[datin$event==TRUE])) # event times
  studyt<-unique(max(datin$time))
  if(ts[1]==0) ts<-ts[-1]
  if(ts[length(ts)]!=studyt) ts<-c(ts,studyt)
  
  set2<-datin[datin$time %in% ts,] # create set with only event times
  set2$tstop<-set2$time
  set2<-set2[,!(names(set2) %in% c("time","intlen"))]
  set2<-set2[order(set2$tstop,set2$ID),]
  
  # Set tstop to event time+1 for coxph
  loc<-match(ts,set2$tstop)
  loc1<-loc[-1]
  loc2<-c(loc[-c(1,2)]-1,length(set2$tstop))
  rm(loc)
  repl<-loc2-loc1+1
  tnew<-rep(ts[-length(ts)], times=repl)
  repz<-as.numeric(length(set2$tstop)-length(tnew))
  set2$time<-c(rep(0,repz),tnew)
  rm(tnew,repz,ts,repl,loc1,loc2)
  
  set2<-set2[!is.na(set2$tstop),]
  return(set2)
}
