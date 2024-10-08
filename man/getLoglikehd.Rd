% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGeneric.R, R/AllMethod.R
\name{getLoglikehd}
\alias{getLoglikehd}
\alias{getLoglikehd,HMMVB-method}
\alias{getLoglikehd,HMMVBBIC-method}
\alias{getLoglikehd,HMMVBclust-method}
\title{Accessor for 'Loglikehd' slot.}
\usage{
getLoglikehd(object)

\S4method{getLoglikehd}{HMMVB}(object)

\S4method{getLoglikehd}{HMMVBBIC}(object)

\S4method{getLoglikehd}{HMMVBclust}(object)
}
\arguments{
\item{object}{Object of class "HMMVB", "HMMVBBIC" "HMMVBclust".}
}
\description{
This function outputs Loglikelihood for each data point in a trained HMM-VB model or Loglikelihood for a new dataset in a HMM-VB model.
}
