# Default Cp.vec, ts.vec, tf.vec, intlen, studyt, and structure already in "makeDVecs" function
ScenSpec=function(half=90,st.dose=1,baser=0.1,OR=1.5,struct=0,Cp.vec=c(1),ts.vec=c(0),tf.vec=c(900),intlen=1,studyt=NULL) {
  
  
  if(is.null(studyt)) studyt=max(tf.vec) #set study length equal to total follow-up time if not explicitly defined
  #Pull in dataframe with time intervals and dosings
  made=makeDVecs(struct=struct,Cp.vec=Cp.vec,ts.vec=ts.vec,tf.vec=tf.vec,intlen=intlen,studyt=studyt)
  
  #For risk assessment, set prevalence rate beta per interval per study length
  frac=intlen/studyt 
  beta0=log(baser*frac)
  #Increased odds of risk as a function of the standard dose of risk
  beta1=log(OR)/st.dose
  
  currC=C1fun.h(thalf=half,dat=made) 
  prob=exp(beta0+beta1*currC)/(1+exp(beta0+beta1*currC))
  
  temp=cbind(half,OR,baser,st.dose,currC,prob,intlen,made)
  temp
}

#Default Cp.vec, ts.vec, tf.vec, intlen, studyt, and structure already in "makeDVecs" function
#Similar to ScenSpec(), but here you can include different half-life specs for incline vs. decline
ScenSpec2=function(half=c(50,100),st.dose=1,baser=0.1,OR=1.5,struct=0,Cp.vec=c(1),ts.vec=c(0),tf.vec=c(900),intlen=1,studyt=NULL) {
  
  if(length(half)==1) {
    print("Assuming the same half-life for incline and decline")
    half<-c(half,half)
  }
  if(length(half)>2) return(print("ERROR: Only one half per direction allowed"))
  
  if(is.null(studyt)) studyt=max(tf.vec) #set study length equal to total follow-up time if not explicitly defined
  #Pull in dataframe with time intervals and dosings
  made=makeDVecs(struct=struct,Cp.vec=Cp.vec,ts.vec=ts.vec,tf.vec=tf.vec,intlen=intlen,studyt=studyt)
  
  #For risk assessment, set prevalence rate beta per interval per study length
  frac=intlen/studyt 
  beta0=log(baser*frac)
  #Increased odds of risk as a function of the standard dose of risk
  beta1=log(OR)/st.dose
  
  currC=C1fun.2h(thalf=half,dat=made)$Conc 
  prob=exp(beta0+beta1*currC)/(1+exp(beta0+beta1*currC))
  
  temp=cbind(half.in=half[1],half.out=half[2],OR,baser,st.dose,currC,prob,intlen,made)
  temp
}
