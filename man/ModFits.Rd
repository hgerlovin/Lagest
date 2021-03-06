% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ModFits.R
\name{fit.cph}
\alias{fit.cph}
\alias{fit.pool}
\title{Fit Cox Proportional Hazards Regression or Pooled Logistic Regression Models}
\usage{
fit.cph(C.in, dat, covs = NULL)

fit.pool(C.in, dat, covs = NULL)
}
\arguments{
\item{C.in}{Vector of values for the effective exposure that corresponds in length and order with the either Surv(time,tstop,event) outcome matrix (fit.cph) or length of the event column in the dataframe \code{dat}.}

\item{dat}{Dataframe with appropriately named columns to be used in either "Surv(time,tstop,event)~" or "event~" and that contains all covariate columns listed in the \code{covs} input. For the fit.cph(), it is assumed that the dataset is already in the appropriate survival dataset format.}

\item{covs}{Vector of names of the covariates to use in the model. Default is none. Function collapses the vector into a formula object separated by "+". For fit.cph() - To specify stratifying or cluster parameters, just include the names in the vector as strata() or cluster() - i.e. covs=c("bmi","strata(age)","hypertension"). Within the find.1half() and find.2half() algorithms, the appropriate strata should be converted to work with this input.}
}
\value{
Outputs a model object corresponding to the function used.
}
\description{
Basic models fit using survival and speedglm packages. For the survival models, Breslow handling of ties is used.
}
