% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/objects.R
\name{processing}
\alias{processing}
\alias{processing,giotto-method}
\title{Preprocess spatial data}
\usage{
processing(object, ...)

\S4method{processing}{giotto}(
  object,
  expression.values = "raw",
  norm.method = c("standard", "osmFISH"),
  min.cells = 10,
  min.features = 200,
  dims = 1:10,
  resolution = 0.6,
  save.plot = FALSE,
  batch.columns = NULL,
  verbose = FALSE,
  outputFolder = NULL,
  python_Path = NULL,
  ...
)
}
\arguments{
\item{object}{An Giotto object}

\item{...}{Argument passed to other methods}

\item{expression.values}{expression values to use}

\item{norm.method}{normalization method to use.}

\item{min.cells}{minimum of cells to expressed in a gene}

\item{min.features}{minimum of genes to detected in a cell}

\item{dims}{number of dimensions to use as input for UMAP}

\item{resolution}{resolution for leiden cluster}

\item{save.plot}{directly save the plot [boolean]}

\item{batch.columns}{metadata columns that represent different batch (max = 2) \code{\link[Giotto]{adjustGiottoMatrix}}}

\item{verbose}{verbose}

\item{outputFolder}{Output folder to save results.}

\item{python_Path}{python path, default as NULL.}
}
\description{
Preprocess spatial data
}
