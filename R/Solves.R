#' Time to reduction in excess hazard
#'
#' \code{solve.time} Calculates the time to \% reduction in the excess hazard following discontinuation of an exposure. \code{solve.red} Determines the time to change from a starting hazard ratio to an ending hazard ratio, given half-lives and plateau of effective exposure.
#' @param beta Effect size coefficient (beta) from the model fitting results.
#' @param D Exposure level prior to discontinuation
#' @param half.in Incline half-life parameter if the resulting model selected was the two-parameter effective exposure model. For the one-parameter effective exposure model, use the same value for the \code{half.in} and \code{half.out} inputs.
#' @param half.out Decline half-life parameter if the resulting model selected was the two-parameter effective exposure model. See \code{half.in} for how to handle OPEE framework.
#' @param time.in Amount of time exposed prior to discontinuation. Value should exceed 1.
#' @param reduct Desired excess hazard reduction. For example, to determine the time to 50\% reduction in the hazard ratio for a specific individual use \code{reduct=0.5}.
#' @param s.val Starting hazard ratio value for \code{solve.red()} function.
#' @param e.val End-point hazard ratio value for \code{solve.red()} function.
#' 
#' @return Returns values for the estimated starting (time=0 at discontinuation) hazard ratio (HR), starting risk or log(HR), ending HR and log(HR), time required to return to the reduced HR, and the relative proportion of reduction from start to end time on both HR and log(HR) scales.
#' \describe{
#'   \item{start.risk}{The starting log(HR) after \code{time.in} units-time of exposure prior to discontinuation.}
#'   \item{start.relrisk}{The starting HR after \code{time.in} units-time of exposure prior to discontinuation.}
#'   \item{time.needed}{The calculated time needed to reduce the excess hazard to \code{reduct}.}
#'   \item{end.risk}{The ending log(HR) for the individual following the \code{time.needed} units-time.}
#'   \item{end.relrisk}{The ending HR for the individual following the \code{time.needed} units-time.}
#'   \item{relrisk.red}{Proportion of reduction in excess hazard. Note that this returns the input \code{reduct} value.}
#'   \item{risk.red}{Proportion of reduction in the log(HR) scale. }
#'   \item{startEE}{Effective Exposure starting value based on the dosing scale with maximum HR at the value of \code{D} input parameter.}
#'   \item{endEE}{Effective Exposure ending value that corresponds to the reduced excess hazard. Similarly, this is relative to the value of \code{D} input parameter.}
#' }
#' @export
#' @examples
#' 
#' ## Calculate the time to 50% reduction in the HR for a 2 packs/day smoker of 30-years. 
#' # Final model being used comes from results in Chapter 4 (OPEE Packs/Day Dosing in Full BWHS Sample).
#' 
#' solve.time(beta=log(2.63), D=2, half.in=5.85, half.out=5.85, time.in=30, reduct=0.5)
#' 
#' ## Calculate the time to 50% reduction in the HR for a 1 pack/day smoker of 30-years. 
#' # Final model being used comes from results in Chapter 4 (OPEE Packs/Day Dosing in Full BWHS Sample).
#' 
#' solve.time(beta=log(2.63), D=1, half.in=5.85, half.out=5.85, time.in=30, reduct=0.5)
solve.time<-function(beta,D,half.in,half.out,time.in,reduct){
  sEE<-D*(1-exp(-log(2)*time.in/half.in)) 
  s.val1<-beta*sEE# determine the Risk starting point
  s.val<-exp(s.val1) # starting relative risk
  e.val<-s.val-((s.val-1)*reduct) # determine percent reduction of the excess risk
  e.val1<-log(e.val) # end risk
  endEE<-e.val1/beta
  risk.red<-abs(s.val1-e.val1)/s.val1
  ratiov<-e.val1/s.val1 # determine the ratio on the EE scale
  lnr<-log(ratiov) # prepare to solve for t
  t<-lnr/(-log(2)/half.out) # divide by the out parameter
  return(c(start.risk=s.val1,start.relrisk=s.val,time.needed=t,end.risk=e.val1,end.relrisk=e.val,relrisk.red=reduct,risk.red=risk.red,startEE=sEE,endEE=endEE))
}


#' @rdname solve.time  
#' @export
#' @examples
#' 
#' ## Calculate the time to change from HR of 6.5 to 2.5, with given OPEE half-life of 5.85 years
#' # for an effective exposure plateau HR of 2.6
#' solve.red(beta=log(2.63),s.val=6.5,e.val=2.6,half.in=5.85,half.out=5.85)
solve.red<-function(beta,s.val,e.val,half.in,half.out){
  s.val1<-log(s.val) # determine the Risk starting point
  sEE<-s.val1/beta
  e.val1<-log(e.val) # end risk
  endEE<-e.val1/beta
  reduct<-abs(s.val-e.val)/(s.val-1) # determine percent reduction of the excess risk
  risk.red<-abs(s.val1-e.val1)/s.val1
  ratiov<-e.val1/s.val1 # determine the ratio on the EE scale
  lnr<-log(ratiov) # prepare to solve for t
  t<-lnr/(-log(2)/half.out) # divide by the out parameter
  return(c(start.risk=s.val1,start.relrisk=s.val,time.needed=t,end.risk=e.val1,end.relrisk=e.val,relrisk.red=reduct,risk.red=risk.red,startEE=sEE,endEE=endEE))
}