#' Simulate single-scenario dataset
#'
#' Generate simulated datasets based on single scenario specification. Can use either one-parameter or two-parameter effective exposure models.
#' @inheritParams ScenSpec 
#' @param scendat Pre-defined scenario to use for simulating subjects (only in DatScen() function). Default assumes that a Scen1 dataframe exists and should be used.
#' @param repper Total number of subjects to simulate. Default is 10,000 subjects.
#' @param seedno Random number generator seed specification for simulating event vs. probability matrix. Default is 25.
#' @return Dataframe with subject-time-incremented observations/rows. 
#' @return Columns taken directly from the pre-simulated scenario trajectory include the input parameter values (st.dose, intlen), values for time, currD, everD, and the true effective exposure (currC). 
#' @return Additional columns created for event, subject ID, and scenario number. Three columns per regimen reflect the exposure-specific dose, time since start, and time since end for each subject-time-specific observation: DoseX, tStartX, tEndX. 
#' @return Each subject will only have one row per time unit up until the first event - iterative comparison of random draw (Unif(0,1)) to conditional logistic probability starts at time 0. 
#' @return Specifics regarding the output dataframe columns (colNames):
#' \describe{
#'    \item{Dose1...DoseX}{Columns indicating the overall doses for each regimen. Repeated throughout for computational ease.}
#'    \item{tStart1...tStartX}{Columns indicating the time since starting the specific regimen -- depends on the point in the trajectory. i.e. Takes a value of 0 for times prior to initiation and increments parallel with time following initiation.}
#'    \item{tEnd1...tEndX}{Columns indicating the time since discontinuing the specific regimen -- depends on the point in the trajectory. i.e. Takes a value of 0 for times prior to start of regimen and while regimen is "on". Increments parallel to time following discontinuation.}
#'    \item{time}{Column for the study time at the observation.}
#'    \item{currD}{Column with value for the current regimen dose for the subject-time-specific observation.}
#'    \item{everD}{Column indicating whether any exposure has occurred as of (prior to and including) the subject-time-specific observation.}
#'    \item{currC}{Column with "true" effective exposure at the given time.}
#'    \item{st.dose}{Input standard dose value, repeated down the column for all time-points.}
#'    \item{intlen}{Input time increment value, repeated down the column for all time-points.}
#'    \item{event}{Column taking value of 1 for subject-time-specific events and 0 otherwise.}
#'    \item{ID}{Subject identifier. For DatScen() -Increments from 1 to input value "repper". Dataframe sorted by ID and then time. For MultiScen() - Increments from 1 to total(nperscen). If one value is used for input to "nperscen", such as 5000, the IDs for subjects in scenario 1 will range from 1 to 5000, scenario 2 from 5001 to 10000, and scenario 3 from 10001 to 15000. For nperscen=c(5000,5000,2000), the IDs for scenario 3 would range from 10001 to 12000. Dataframe sorted by ID and then time.}
#'    \item{scen.num}{Scenario number. Not particularly useful in output from DatScen(). For MultiScen() only - Subjects generated in the order of the scenario specification -- this means that changing the order in "combo" would alter the RNG matrix comparisons.}
#' }
#' @describeIn DatScen Assumes a pre-generated scenario dataframe that feeds into the "scendat" value. Underlying dataframe can come from either OPEE- or TPEE-generated scenarios.
#' @seealso ScenSpec, ScenSpec2, MultiScen
#' @export
DatScen <- function(scendat=Scen1, repper=10000, seedno=25) {
  time=scendat[,"time"] #Check .subset2 function to speed things up #looked it up and could not find a faster alternative
  prob=scendat[,"prob"]
  
  keeps=scendat[,!(names(scendat) %in% c("half","prob","OR","baser"))] #requires the names(scendat)
  
  set.seed(seedno)
  
  event=matrix(NA,nrow=repper,ncol=length(time))
  r=matrix(runif(length(time)*repper),nrow=repper,ncol=length(time))
  for(i in 1:repper) {
    for(t in 1:length(time)){
      event[i,t]=(r[i,t]<prob[t])
      if(event[i,t]==TRUE)  break
    }
  }
  newdat<-reshape::melt.matrix(event)
  colnames(newdat)<-c("ID","time","event")
  newdat=newdat[!is.na(newdat$event),]
  newdat=merge(newdat,keeps,by=c("time"))
  
  newdat=newdat[with(newdat, order(ID,time)),]
  rownames(newdat)=NULL
  newdat$scen.num<-1
  
  return(newdat)
}

#' @describeIn DatScen To input the OPEE scenario specifics from which to generate a random sample. Uses ScenSpec() to generate the underlying scenario. Not recommended - please pre-generate scenarios and use either DatScen or MultiScen functions instead.
#' @export
DatSpec <- function(repper=10000,half=90,st.dose=1,baser=0.1,OR=1.5,struct=0,Cp.vec=c(1),ts.vec=c(0),tf.vec=c(900),intlen=1,studyt=NULL,seedno=25) {
  
  Scen1<-ScenSpec(half=half,st.dose=st.dose,baser=baser,OR=OR,struct=struct,Cp.vec=Cp.vec,ts.vec=ts.vec,tf.vec=tf.vec,intlen=intlen,studyt=studyt)
  return(DatScen(scendat=Scen1,repper=repper,seedno=seedno))
}

#' @describeIn DatScen To input the TPEE scenario specifics from which to generate a random sample. Uses ScenSpec2() to generate the underlying scenario. Not recommended - please pre-generate scenarios and use either DatScen or MultiScen functions instead.
#' @export
DatSpec2 <- function(repper=10000,half=c(50,100),st.dose=1,baser=0.1,OR=1.5,struct=0,Cp.vec=c(1),ts.vec=c(0),tf.vec=c(900),intlen=1,studyt=NULL,seedno=25) {
  
  Scen1<-ScenSpec2(half=half,st.dose=st.dose,baser=baser,OR=OR,struct=struct,Cp.vec=Cp.vec,ts.vec=ts.vec,tf.vec=tf.vec,intlen=intlen,studyt=studyt)
  return(DatScen(scendat=Scen1,repper=repper,seedno=seedno))
}
