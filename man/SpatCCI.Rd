% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cellphonedb.r
\name{SpatCellphoneDB}
\alias{SpatCellphoneDB}
\alias{cellphonedb_merge}
\title{cellphoneDB for spot-population}
\usage{
SpatCellphoneDB(
  object,
  outputFolder = NULL,
  numCores = 4,
  Geneformat = "hgnc_symbol",
  cellphonedbPath = NULL,
  method = "kmeans",
  ...
)

cellphonedb_merge(interaction)
}
\arguments{
\item{object}{A giotto object}

\item{outputFolder}{Output folder to save results.}

\item{numCores}{The number of cores to use, default as 4.}

\item{Geneformat}{Format of gene name in count matrix, one of "ensembl_gene", "gene_symbol".}

\item{cellphonedbPath}{Path to launch cellphoneDB.}

\item{method}{Spatial clustering method to annotate spots population, either "kmeans" or "HMRF".}

\item{...}{other arguments}

\item{interaction}{A list out of cellphoneDB calculation.}
}
\description{
cellphoneDB for spot-population

cellphonedb_merge
}
