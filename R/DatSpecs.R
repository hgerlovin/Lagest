DatSpec <- function(repper=10000,half=90,st.dose=1,baser=0.1,OR=1.5,struct=0,Cp.vec=c(1),ts.vec=c(0),tf.vec=c(900),intlen=1,studyt=NULL,seedno=25) {
  
  Scen1=ScenSpec(half=half,st.dose=st.dose,baser=baser,OR=OR,struct=struct,Cp.vec=Cp.vec,ts.vec=ts.vec,tf.vec=tf.vec,intlen=intlen,studyt=studyt)
  return(DatScen(scendat=Scen1,repper=repper,seedno=seedno))
}

DatSpec2 <- function(repper=10000,half=c(50,100),st.dose=1,baser=0.1,OR=1.5,struct=0,Cp.vec=c(1),ts.vec=c(0),tf.vec=c(900),intlen=1,studyt=NULL,seedno=25) {
  
  Scen1=ScenSpec2(half=half,st.dose=st.dose,baser=baser,OR=OR,struct=struct,Cp.vec=Cp.vec,ts.vec=ts.vec,tf.vec=tf.vec,intlen=intlen,studyt=studyt)
  return(DatScen(scendat=Scen1,repper=repper,seedno=seedno))
}

DatScen <- function(scendat=Scen1, repper=10000, seedno=25) {
  require(reshape)
  time=scendat[,"time"] #Check .subset2 function to speed things up #looked it up and could not find a faster alternative
  prob=scendat[,"prob"]
  
  keeps=scendat[,!(names(scendat) %in% c("half","prob","OR","baser"))] #!(names %in% c("half","prob","currC")) #requires the names(scendat)
  
  set.seed(seedno)
  
  event=matrix(NA,nrow=repper,ncol=length(time))
  r=matrix(runif(length(time)*repper),nrow=repper,ncol=length(time))
  for(i in 1:repper) {
    for(t in 1:length(time)){
      event[i,t]=(r[i,t]<prob[t])
      if(event[i,t]==TRUE)  break
    }
  }
  newdat=melt.matrix(event)
  colnames(newdat)<-c("ID","time","event")
  newdat=newdat[!is.na(newdat$event),]
  newdat=merge(newdat,keeps,by=c("time"))
  
  newdat=newdat[with(newdat, order(ID,time)),]
  rownames(newdat)=NULL
  newdat$scen.num<-1
  
  return(newdat)
}
