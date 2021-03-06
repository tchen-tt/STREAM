% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/objects.R
\name{CCI}
\alias{CCI}
\alias{CCI,giotto-method}
\title{Cell and cell interaction}
\usage{
CCI(object, ...)

\S4method{CCI}{giotto}(object, fdrCutoff = 0.1, numGens = 20, numCores = 8, ...)
}
\arguments{
\item{object}{a giotto object}

\item{...}{other parameters for scHOT function \code{\link[scHOT]{scHOT}}}

\item{fdrCutoff}{FDR value cutoff}

\item{numGens}{number of genes for gene paire}

\item{numCores}{The number of cores to use}
}
\description{
Cell and cell interaction

gene pair coexpression based on scHOT
}
