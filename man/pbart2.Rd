% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bart_functions.R
\name{pbart2}
\alias{pbart2}
\title{Convenience wrapper for pbart from the BART package.}
\usage{
pbart2(x.train, y.train, x.test = NULL, mc.cores = 1, ndpost = 1000,
  keeptrainfits = FALSE, verbose = FALSE, ...)
}
\arguments{
\item{x.train}{Data frame containing the X variables for the training data}

\item{y.train}{Vector containing the Y variable for the training data}

\item{x.test}{Xata frame containing the X variables for the test data}

\item{mc.cores}{Number of cores to use}

\item{ndpost}{Number of posterior samples to draw}

\item{keeptrainfits}{Keep samples for training data?}

\item{verbose}{Show pbart progress messages?}

\item{...}{Additional BART parameters passed to \code{pbart}}
}
\value{
Returns a pbart object
}
\description{
Convenience wrapper for pbart from the BART package.
}
