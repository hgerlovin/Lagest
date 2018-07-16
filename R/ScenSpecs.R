#' OPEE/TPEE dataset creation
#'
#' Generates an OPEE trajectory based on pre-specified inputs. Output dataframe can be used in DatScen and MultiScen.
#' @param half One-parameter half-life for simulating the relative concentration/effective exposure, currC. Default assumes 90 time unit half-life.
#' @param st.dose Relative units for dose plateau. Default assumes binary exposure.
#' @param baser Overall prevalence input. Parameter feeds into the intercept beta via \eqn{log(baser\times(intlen/studyt))}. This beta is then used to calculate the subject-time-specific conditional logistic probability = \eqn{\dfrac{exp(beta0 + beta1*currC)}{(1 + exp(beta0 + beta1*currC))}}. Default is 10% prevalence of the outcome throughout the course of the study follow-up.
#' @param OR Proportional Odds/Hazard Ratio. Parameter feeds into the beta1 via ln(OR/st.dose). This beta is then used to calculate the subject-time-specific conditional logistic probability = \eqn{\dfrac{exp(beta0 + beta1*currC)}{(1 + exp(beta0 + beta1*currC))}}. Default assumes that steady state full exposure increases the odds of event 50%.
#' @param struct Structure indicator. If turned on (1), then additional regimen is added for time following discontinuation. Default is off (0), assuming that the total number of regimens is fixed and does not need additional follow-up.
#' @param Cp.vec Vector of doses for the regimens. Default assumes the binary exposure plateau and that there is a single regimen of exposure (i.e. \code{Cp.vec=c(1)}). To have multiple exposures, include the same number of vector components in \code{Cp.vec}, \code{ts.vec}, and \code{tf.vec}.
#' @param ts.vec Vector of start times for the regimens. Default assumes the exposure was started at time 0.
#' @param tf.vec Vector of end times for the regimens. Default assumes the exposure continues through time=900. When \code{studyt} is not specified, the last specified end-time (last regimen) is used as the total study time.
#' @param intlen Increment time to use. Default is 1 time unit.
#' @param studyt Total study follow-up time. Default is \code{NULL} and will pull the last regimen stop time.
#' 
#' @return Outputs dataframe of study intervals with corresponding/useable dose regimen vectors, current dose vector.
#' 
#' @param \code{Dose1}...\code{DoseX}	Columns indicating the overall doses for each regimen. Repeated throughout for computational ease. 
#' @param \code{tStart1}...\code{tStartX}	Columns indicating the time since starting the specific regimen -- depends on the point in the trajectory. i.e. Takes a value of 0 for times prior to initiation and increments parallel with time following initiation.
#' @param \code{tEnd1}...\code{tEndX}	Columns indicating the time since discontinuing the specific regimen -- depends on the point in the trajectory. i.e. Takes a value of 0 for times prior to start of regimen and while regimen is "on". Increments parallel to time following discontinuation.
#' @param \code{time}	Column for the study time at the observation.
#' @param \code{currD}	Column with value for the current regimen dose for the subject-time-specific observation.
#' @param \code{everD}	Column indicating whether any exposure has occurred as of (prior to and including) the subject-time-specific observation.
#' @param \code{currC}	Column with "true" effective exposure at the given time.
#' @param \code{prob}	Column with time-specific conditional logistic probability of event. This is the value used for assigning events during data generation - compared to random draw from Unif(0,1).
#' @param \code{half}	Input half-life value, repeated down the column for all time-points. Retained for simulation purposes and later discarded.
#' @param \code{OR}	Input Odds Ratio value, repeated down the column for all time-points. Retained for simulation purposes and later discarded.
#' @param \code{baser}	Input baseline prevalence value, repeated down the column for all time-points. Retained for simulation purposes and later discarded.
#' @param \code{st.dose}	Input standard dose value, repeated down the column for all time-points. Retained for simulation purposes and later discarded.
#' @param \code{intlen}	Input time increment value, repeated down the column for all time-points. Retained for simulation purposes and later discarded.
#' 
#' @export
#' 
#' @example 
#' ScenSpec()


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

#' @describeIn ScenSpec Generates an TPEE trajectory based on pre-specified inputs.
#' @export

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
