% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HessDerivFun2.R
\name{timesums1}
\alias{timesums1}
\alias{timesums2}
\title{Hessian-derivation functions}
\usage{
timesums1(b1, Ct, C2, C3, y, p)

timesums2(b1, Ct, dC.1, dC.2, d2C.1, d2C.2, d2C.12, y, p)
}
\arguments{
\item{b1}{Final beta parameter estimate using the estimated half-life (or half-lives)}

\item{Ct}{Vector or column of dataframe with the calculated effective exposure values}

\item{C2, dC.1, dC.2}{Vector or column of dataframe with the calculated first derivative of the effective exposure. C2 is used in the OPEE framework, while dC.1 and dC.2 correspond to the first derivatives of the TPEE model framework effective exposure with respect to the incline and decline half-lives, respectively.}

\item{C3}{Vector or column of dataframe with the calculated second derivative of the effective exposure. C3 is only applicable in the OPEE framework.}

\item{y}{A vector or column in the dataframe with the true event values. Should be the same length and correspond to the observations for the other input parameters. Takes a value of 1 for events and 0 for non-events.}

\item{p}{The predicted probability of event for a given subject-time-observation using the fully-specified model.}

\item{d2C.1, d2C.2, d2C.12}{Vector or columns of the dataframe with the calculated second derivatives of the effective exposure values with respect to the incline, decline, and both parameters, respectively. Applicable to the TPEE framework only.}
}
\value{
Output includes:
\describe{
   \item{u}{Score for the model fit under given parameters}
   \item{I}{Fisher's Information for the model fit under given parameters}
   \item{llik}{Log-likelihood for the model fit under given parameters}
}
}
\description{
Functions here are for calculating the first and second derivatives for the 1 and 2 lag parameter Effective Exposure. Used in the half-life finding algorithms to compute SE. Continuation from the initial Hessian-derivation functions.
}
\section{Functions}{
\itemize{
\item \code{timesums1}: Calculate the Hessian matrix components from the subject-time-specific equation values for the CPH OPEE model

\item \code{timesums2}: Calculate the Hessian matrix components from the subject-time-specific equation values for the CPH TPEE model
}}

