% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MultiScen.R
\name{MultiScen}
\alias{MultiScen}
\alias{makeDat.cph}
\title{Simulate multiple-scenarios dataset}
\usage{
MultiScen(combo = list(scen1, scen2, scen3), nperscen = 5000, seedno = 25)

makeDat.cph(combo = list(base1.up, base1.down, base1.ctrl), nperscen = 5000,
  seedno = seedset)
}
\arguments{
\item{combo}{Pre-defined scenarios to use for simulating subjects. Must be input as a list of ScenSpec/ScenSpec2 output objects.}

\item{nperscen}{Total number of subjects to simulate per trajectory/scenario. Default is 5,000 subjects. Can be input as a vector of values with the same number of elements as listed in "combo" - i.e. nperscen=c(5000,5000,2000) would imply 5,000 subjects generated under scenarios 1 and 2 and only 2,000 subjects with scenario 3.}

\item{seedno}{Random number generator seed specification for simulating event vs. probability matrix. Default is 25.}
}
\value{
Dataframe with subject-time-incremented observations/rows.

Columns taken directly from the pre-simulated scenario trajectory include the input parameter values (st.dose, intlen), values for time, currD, everD, and the true effective exposure (currC).

Additional columns created for event, subject ID, and scenario number. Three columns per regimen reflect the exposure-specific dose, time since start, and time since end for each subject-time-specific observation: DoseX, tStartX, tEndX.

Each subject will only have one row per time unit up until the first event - iterative comparison of random draw (Unif(0,1)) to conditional logistic probability starts at time 0.

Specifics regarding the output dataframe columns (colNames):
\describe{
   \item{Dose1...DoseX}{Columns indicating the overall doses for each regimen. Repeated throughout for computational ease.}
   \item{tStart1...tStartX}{Columns indicating the time since starting the specific regimen -- depends on the point in the trajectory. i.e. Takes a value of 0 for times prior to initiation and increments parallel with time following initiation.}
   \item{tEnd1...tEndX}{Columns indicating the time since discontinuing the specific regimen -- depends on the point in the trajectory. i.e. Takes a value of 0 for times prior to start of regimen and while regimen is "on". Increments parallel to time following discontinuation.}
   \item{time}{Column for the study time at the observation.}
   \item{currD}{Column with value for the current regimen dose for the subject-time-specific observation.}
   \item{everD}{Column indicating whether any exposure has occurred as of (prior to and including) the subject-time-specific observation.}
   \item{currC}{Column with "true" effective exposure at the given time.}
   \item{st.dose}{Input standard dose value, repeated down the column for all time-points.}
   \item{intlen}{Input time increment value, repeated down the column for all time-points.}
   \item{event}{Column taking value of 1 for subject-time-specific events and 0 otherwise.}
   \item{ID}{Subject identifier. For DatScen() -Increments from 1 to input value "repper". Dataframe sorted by ID and then time. For MultiScen() - Increments from 1 to total(nperscen). If one value is used for input to "nperscen", such as 5000, the IDs for subjects in scenario 1 will range from 1 to 5000, scenario 2 from 5001 to 10000, and scenario 3 from 10001 to 15000. For nperscen=c(5000,5000,2000), the IDs for scenario 3 would range from 10001 to 12000. Dataframe sorted by ID and then time.}
   \item{scen.num}{Scenario number. Not particularly useful in output from DatScen(). For MultiScen() only - Subjects generated in the order of the scenario specification -- this means that changing the order in "combo" would alter the RNG matrix comparisons.}
}
}
\description{
Generate simulated datasets based on multiple pre-specified scenarios. Should only use scenarios generated using the same underlying effective exposure structure (same number of parameters and values).
}
\section{Functions}{
\itemize{
\item \code{makeDat.cph}: Used for parallel processing in simulations to make a survival-based dataset following the MultiScen() function. This step removes all time-point observations that do not have a corresponding event - i.e. Only unique event-times are retained. This combines the MultiScen() and full.surv() functions to output a single dataset for simulation purposes.
}}

\seealso{
ScenSpec, ScenSpec2, DatSpec, DatSpec2, DatScen
}
