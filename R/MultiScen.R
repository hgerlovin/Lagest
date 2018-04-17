# Reads in pre-created scenarios, set number of subjects per scenario, set random seed
# Creates dataset for nperscen subjects for each scenario in the combo list using given seed
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
    
    newerdat[[snum]]<-melt.matrix(events[[snum]])
    
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
  
  fin.dat<-do.call(rbind.fill,newerdat)
  rownames(fin.dat)=NULL
  
  return(fin.dat)
}