#' Simplified version of constructing dose regimen vectors.
#'
#' Reads in dose, start and end times, and half-life. Calculates concentration for inputs.
#' @param Cp.vec Vector of doses for the regimens. Default assumes the binary exposure plateau and that there is a single regimen of exposure. To have multiple exposures, include the same number of vector components in \code{Cp.vec}, \code{ts.vec}, and \code{tf.vec}.
#' @param ts.vec Vector of start times for the regimens. Default assumes the exposure was started at time 0.
#' @param tf.vec Vector of end times for the regimens. Default assumes the exposure continues through time=900. When \code{studyt} is not specified, the last specified end-time (last regimen) is used as the total study time.
#' @param intlen Increment time to use. Default is 1 time unit.
#' @param studyt Total study follow-up time. Default is \code{NULL} and will pull the last regimen stop time.
#' @param structure Structure indicator. If turned on (1), then additional regimen is added for time following discontinuation. Default is off (0), assuming that the total number of regimens is fixed and does not need additional follow-up.
#' @return Outputs dataframe of study intervals with corresponding/useable dose regimen vectors, current dose vector.
#' 
#' \code{Dose1}...\code{DoseX} Columns indicating the overall doses for each regimen. Repeated throughout for computational ease.
#' 
#' \code{tStart1}...\code{tStartX} Columns indicating the time since starting the specific regimen -- depends on the point in the trajectory. i.e. Takes a value of 0 for times prior to initiation and increments parallel with time following initiation.
#' 
#' \code{tEnd1}...\code{tEndX}  Columns indicating the time since discontinuing the specific regimen -- depends on the point in the trajectory. i.e. Takes a value of 0 for times prior to start of regimen and while regimen is "on". Increments parallel to time following discontinuation.
#' 
#' \code{time} & Column for the study time at the observation.
#' 
#' \code{currD} & Column with value for the current regimen dose for the subject-time-specific observation.
#' 
#' \code{everD} &	Column indicating whether any exposure has occurred as of (prior to and including) the subject-time-specific observation.
#' @export
#' @examples
#' makeDVecs()


makeDVecs=function(struct=0,Cp.vec=c(1),ts.vec=c(0),tf.vec=c(900),intlen=1,studyt=NULL) {
  
  nCp=length(Cp.vec) #Set number of dosing levels
  
  #Pre-specifying that study end time is dependent on decay from dosing. Only use this for discontinued drugs with follow-up
  if(struct==1 & Cp.vec[nCp]!=0) { #If final dose is not equal to 0, set to 0 with follow-up time for 20 half-lives
    nCp=nCp+1
    Cp.vec[nCp]=0
    tf.vec[nCp]=half*20
    ts.vec[nCp]=tf.vec[nCp-1]
  }
  
  if(is.null(studyt)) studyt=max(tf.vec) #set study length equal to total follow-up time if not explicitly defined
  time=seq(0,to=studyt,by=intlen)
  
  #Create lists of values for each individual
  delTstart=lapply(ts.vec, function(x) unlist(lapply(time, function(t) max(0,t-x)),use.names=F))
  delTend=lapply(tf.vec, function(x) unlist(lapply(time, function(t) max(0,t-x)),use.names=F))
  Dose=lapply(Cp.vec,function(q) rep(q,length=length(time)))
  
  names(delTend)=paste0("tEnd",1:nCp)
  names(delTstart)=paste0("tStart",1:nCp)
  names(Dose)=paste0("Dose",1:nCp)
  
  currD=colSums(do.call(rbind,lapply(1:nCp, function(q) (time<=tf.vec[q] & time>ts.vec[q])*Cp.vec[q])))
  everD=colSums(do.call(rbind,lapply(1:nCp, function(t) ifelse(any(Cp.vec[1:t]>0),1,0))))
  
  temp=as.data.frame(cbind(time,currD,everD,as.data.frame(Dose),as.data.frame(delTstart),as.data.frame(delTend)))
  temp=temp[(time<=studyt),] #For study lengths shorter than dosing length, remove extra measurement intervals
  temp
}
