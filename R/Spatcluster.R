#' @title Kmeans methods to calculate spatial gene expression pattern
#' @param object A giotto object
#' @param cellType Column in meta to assign cell-type-annotation for spatial clustering. Default as "leiden" when external annotation is missing.
#' @param method The distance measure to be used. This must be one 
#' of "euclidean", "maximum", "manhattan", "canberra", "binary" or 
#' "minkowski". Any unambiguous substring can be given.
#' @param k Defines k for the k-nearest neighbor algorithm
#' @param centers The number of clusters
#' @param iter.max The maximum number of iterations allowed
#' @param nstart If centers is a number, how many random sets should be chosen
#' @param ... other arguments \code{\link[stats]{kmeans}}
#'
#' @export
#' @rdname SpatCluster

KmeansCluster <- function(object, cellType = "leiden", method = "maximum",
                          k = 10, centers = 6, iter.max = 300, nstart = 10,
                          ... ) {
  if (!(cellType %in% colnames(object@cell_metadata))) {
    stop("cellType must in cell_metadata", cell. = TRUE)
  }
  method <- match.arg(method, choices = c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"))
  
  if (k > nrow(object@cell_metadata)) {
    k <- nrow(object@cell_metadata) - 1
  }
  
  metaData <- as.data.frame(object@cell_metadata)
  metaData <- metaData[, cellType, drop = TRUE]
  metaData <- cbind.data.frame(metaData, object@spatial_locs)
  # dimnames(metaData) <- list(metaData[, "cell_ID", drop = TRUE], 
  #                           c("cellType", "x", "y", "cell_ID"))
  rownames(metaData) <- metaData[, "cell_ID", drop = TRUE]
  colnames(metaData) <- c("cellType", "x", "y", "cell_ID")
  
  allDist <- dist(metaData[, c("x", "y")], method = method)
  allDist <- as.matrix(allDist)
  
  spotNN <- apply(allDist, 2, function(x) {
    names(x)[order(x)[2:(k+1)]]
  })
  names(spotNN) <- rownames(allDist)
  
  cellType <- unique(metaData[, "cellType", drop = TRUE])
  nCluster <- length(cellType)
  
  propMatrix <- sapply(1:ncol(allDist), function(x) {
    cellTypeCount <- table(metaData[spotNN[, x], "cellType", drop = TRUE])
    allCellCounts <- structure(rep(0L, length = nCluster), names = cellType)
    allCellCounts[match(names(cellTypeCount), cellType)] <- cellTypeCount
    allCellCounts
  })
  propMatrix <- t(propMatrix)
  
  kmeansNF <- kmeans(propMatrix , centers = centers, iter.max = iter.max, nstart = nstart, ...)
  object@cell_metadata$kmeans <- kmeansNF$cluster
  return(object)
}

#' @title HmrfCluster
#' 
#' @param object A giotto object
#' @param expressValue Expression values to use (default [scaled])
#' @param spatialGenes Spatial genes to use for HMRF
#' @param outputFolder Output folder to save results
#' @param kdomain Number of HMRF domains, default 5
#' @param knn number of nearest neighbors based on physical distance
#' @param ... Other arguments. \code{\link{createSpatialKNNnetwork}} \code{\link{doHMRF}} \code{\link{addHMRF}}
#' @rdname SpatCluster
#' @importFrom  Giotto doHMRF addHMRF
#' 
#' @return giotto object
#' @export

HmrfCluster <- function(object, knn = 5, expressValue = "scaled",
                        spatialGenes = NULL, outputFolder = NULL,
                        kdomain = 5,
                        ...) {
  expressValue <- match.arg(expressValue, choices = c("normalized", "scaled", "custom"))
  
  if (is.null(spatialGenes)) stop("spatialGenes should be given.", cell. = TRUE)
  if(is.null(outputFolder)&!is.null(object@instructions$save_dir)){
    outputFolder=object@instructions$save_dir
  }
  if (is.null(outputFolder)) stop("outputFolder should be given.", cell. = TRUE)
  
  object <- createSpatialNetwork(gobject = object, k = knn, name = "D")
  object = createSpatialNetwork(gobject = object, minimum_k = 2, maximum_distance_delaunay = 400)
  
  hmrfSpatialGenes <- doHMRF(gobject = object, expression_values = expressValue, 
                             spatial_genes = spatialGenes, output_folder = outputFolder,
                             k = kdomain, spatial_network_name = "Delaunay_network",overwrite_output = TRUE)
  object <- addHMRF(gobject = object, HMRFoutput = hmrfSpatialGenes, k = kdomain, betas_to_add=50,  ...)
  return(object)
}





