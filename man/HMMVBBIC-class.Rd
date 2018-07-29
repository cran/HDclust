\docType{class}
\name{HMMVBBIC-class}
\alias{HMMVBBIC-class}
\alias{HMMVBBIC}
\alias{show,HMMVBBIC-method}

\title{Class "HMMVBBIC" to represent results of HMM-VB model selection.}
\description{
An S4 class to represent results of HMM-VB model selection. New instances of the class are created by \code{\link{hmmvbBIC}}.
}
\section{Methods}{
\itemize{
\item \bold{show} signature(object = "HMMVBBIC") : show optimal model.
\item \bold{getBIC} signature(object = "HMMVBBIC") : accessor for 'BIC' slot.
\item \bold{getNumst} signature(object = "HMMVBBIC") : accessor for 'numst' slot.
\item \bold{getConfigList} signature(object = "HMMVBBIC") : accessor for 'configList' slot.
\item \bold{getOptHMMVB} signature(object = "HMMVBBIC") : accessor for 'optHMMVB' slot.
}}
\section{Slots}{

\describe{
\item{\code{numst}}{An integer vector specifying the number of mixture components (states) in each variable block for which BIC was calculated during model selection (see \code{\link{hmmvbBIC}}).}

\item{\code{BIC}}{A numeric vector specifying calculated BIC values.}

\item{\code{configList}}{A list of integer vectors specifying the number of states in each variable block for which BIC was calculate (see \code{\link{hmmvbBIC}}).}

\item{\code{optHMMVB}}{The optimal HMM-VB model with smallest BIC value.}
}
}
\seealso{
\code{\link{hmmvbBIC}}
}
