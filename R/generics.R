#' Spatial single cell quality control
#'
#' @param object An Giotto object
#' @param ... Arguments passed to other methods
#'
#' @rdname CellQC
#' @export CellQC
#'

setGeneric(name = "CellQC", def = function(object, ...) object)



#' Preprocess spatial data
#' 
#' 
#' @param object An Giotto object
#' @param ... Argument passed to other methods
#'
#'
#' @rdname processing
#' @export processing

setGeneric(name = "processing", def = function(object, ...) object)


#' Different expression with spatial single
#' 
#' @param object An Giotto object
#' @param ... Argument passed to other methods
#' 
#' @rdname DiffGenes
#' @export
setGeneric(name = "DiffGenes", def = function(object, ...) object)


#' Cluster analysis with spatial single RNA-seq data
#' 
#' @param object An Giotto object
#' @param ... Argument passed to other methods
#' @rdname SpatCluster
#' @export
setGeneric(name = "SpatCluster", def = function(object, ...) object)

#' Cell and cell interaction
#' 
#' @param object An giotto object
#' @param ... Argument passed to other methods
#' 
#' @rdname CCI
#' @export
setGeneric(name = "CCI", def = function(object, ...) object)


#' Integration single cell RNA-seq data and spatial single cell RNA-seq data
#' 
#' @param object An giotto object
#' @param ... Argument passed to other methods
#' 
#' @rdname Integration
#' @export
setGeneric(name = "Integration", def = function(object, ...) object)


#' Deconvolution of Spatial data
#' 
#' @param object An giotto object
#' @param ... Argument passed to other methods
#' 
#' @rdname Deconvolution
#' @export
setGeneric(name = "Deconvolution", def = function(object, ...) object)


#' Visualization of Spatial Clustering
#' 
#' 
#' @param object An Giotto object
#' @param ... Argument passed to other methods
#'
#'
#' @rdname DomainPlot
#' @export DomainPlot

setGeneric(name = "DomainPlot", def = function(object, ...) object)
