% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bayesian_assessment.R
\name{bayes_ms}
\alias{bayes_ms}
\title{Title}
\usage{
bayes_ms(x, y, Pmc = 0.5, optim_details = TRUE, event = 1, ...)
}
\arguments{
\item{x}{a numeric vector of probabilities to be LLO-adjusted. Must only
contain values in [0,1].}

\item{y}{a numeric vector of outcomes corresponding to probabilities in \code{x}.
Must only contain 0 or 1.}

\item{Pmc}{}

\item{optim_details}{Logical.  If \code{TRUE}, the list returned by \code{optim()} when
minimizing the negative log likelihood is also returned by \code{llo_lrt()}.}

\item{event}{Value in \code{y} that represents an "event".  Default value is 1.}

\item{...}{Additional arguments to be passed to \code{optim()}.}
}
\description{
Title
}