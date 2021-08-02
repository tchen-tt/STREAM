#' @include generics.R
#' @importClassesFrom Seurat Seurat
#' @importClassesFrom Giotto giotto
#' @importClassesFrom SPARK SPARK
#' @importClassesFrom scHOT scHOT
#' 
NULL


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Giotto generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' @param expression.values expression values to use
#' @param min.cells 	minimum of cells to expressed in a gene
#' @param min.features minimum of genes to detected in a cell
#' @param norm.method normalization method to use.
#' @param save.plot directly save the plot [boolean]
#' @param batch.columns metadata columns that represent different batch (max = 2) \code{\link[Giotto]{adjustGiottoMatrix}}
#' @param dims number of dimensions to use as input for UMAP
#' @param outputFolder Output folder to save results.
#' @param python_Path python path, default as NULL.
#' @param resolution resolution for leiden cluster
#' @param verbose verbose
#' @importFrom  Giotto installGiottoEnvironment
#' @importFrom Giotto filterGiotto normalizeGiotto addStatistics adjustGiottoMatrix
#' @importFrom Giotto runPCA runUMAP createNearestNetwork doLeidenCluster
#' @import dplyr
#' @rdname processing 
#' @method processing giotto

setMethod("processing", signature = "giotto",
          definition = function(object,
                                expression.values ="raw",
                                norm.method = c("standard", "osmFISH"),
                                min.cells = 10,
                                min.features = 200,
                                dims = 1:10,
                                resolution = 0.6,
                                save.plot = FALSE,
                                batch.columns = NULL,
                                verbose = FALSE,
                                outputFolder=NULL,
                                python_Path=NULL,
                                ...) {
   expression.values <- match.arg(expression.values, choices = c("raw", "normalized", "scaled", "custom"))

   object@spatial_locs$total_counts <- colSums(as.matrix(object@raw_exprs))
   
   
   if(is.null(python_Path)) {
      installGiottoEnvironment()
   }
   
   instrs = createGiottoInstructions(show_plot = FALSE,save_plot = FALSE, 
                                     save_dir = outputFolder,python_path = python_Path)
   
   object@instructions=instrs
                                               
   object <- filterGiotto(gobject = object,
                          expression_values = expression.values,
                          gene_det_in_min_cells = min.cells,
                          min_det_genes_per_cell = min.features,
                          verbose = verbose,
                          ...)
   norm.method <- match.arg(norm.method, choices = c("standard", "osmFISH"))
   object <- normalizeGiotto(gobject = object,
                             norm_methods = norm.method,
                             verbose = verbose)
   
   object <- addStatistics(gobject = object,
                           return_gobject = TRUE)
   
   object <- calculateHVG(gobject = object, 
                          method = "cov_groups",
                          save_plot = save.plot,
                          return_plot = FALSE)
   
   if (!is.null(batch.columns))
      object <- adjustGiottoMatrix(gobject = object)
   
   object <- runPCA(gobject = object, 
                    genes_to_use = Features(object))
   
   object <- runUMAP(gobject = object,
                     dimensions_to_use = dims)
   
   object <- createNearestNetwork(gobject = object,
                                  dimensions_to_use = dims)
   
   object <- doLeidenCluster(gobject = object, 
                             resolution = resolution,
                             name = "leiden")
   
   if (save.plot){
     
     pUMAP <- plotUMAP(gobject = object,cell_color = 'leiden', show_NN_network = F, point_size = 6000/dim(object@spatial_locs)[1], save_plot=FALSE, return_plot=TRUE)
     pUMAP <- pUMAP + scale_colour_npg(alpha=0.4)
     pUMAP <- pUMAP + theme(legend.text = element_text(color="azure4", size = 10, face = "bold"),legend.key.size = unit(0.25, "inches")) +
        guides(colour = guide_legend(override.aes = list(size = 4)))
     ggsave(filename=paste0(outputFolder,"/UMAP.png"),plot = pUMAP,width = 8,height = 5)
     
     markers_scarn = findMarkers_one_vs_all(gobject=object, method="scran", expression_values="normalized", cluster_column="leiden", min_genes=5)
     markertable = markers_scarn %>% select("cluster","genes") %>%
     mutate(num=seq(1,dim(markers_scarn)[1]))%>% 
     group_by(cluster) %>% mutate(row_number = row_number(num))
     markergenes_scran = unique(markertable$genes[which(markertable$row_number<=8)])
     heatmap = plotMetaDataHeatmap(object, expression_values="normalized", metadata_cols=c("leiden"), selected_genes=markergenes_scran, return_plot = TRUE,save_plot = FALSE)
     heatmap = heatmap + 
        theme_classic() + 
        theme(axis.line = element_blank(), 
              axis.ticks = element_blank(), 
              axis.text.y = element_text(face = "bold",size = 4), 
              axis.title.y = element_blank()) + 
        labs(x = "leiden cluster") + 
        coord_equal(ratio=2*length(unique(object@cell_metadata$leiden))/length(markergenes_scran))
     ggsave(filename = paste0(outputFolder,"/scranMarkerHeatmap.png"),plot=heatmap,height = 8)
   }
   object
})


#' @param object A SPARK object
#' @param num.cores The number of cores to use, default 4.
#' @param verbose Output fitting information, default [FALSE]
#' @param ... Other graguments in \code{\link[SPARK]{spark.vc}} 
#' \code{\link[SPARK]{spark.test}}
#' 
#' @return Different expression data.fame
#' @importFrom SPARK spark.vc spark.test
#' @rdname DiffGenes
#' @method DiffGenes SPARK

setMethod("DiffGenes", signature = "SPARK",
          definition = function(object,
                                num.cores = 4,
                                verbose = FALSE,
                                ...) {
   object@lib_size <- Matrix::colSums(object@counts, na.rm = TRUE) # fast fast !!!!
   
   object <- spark.vc(object = object,
                      lib_size = object@lib_size,
                      num_core = num.cores,
                      verbose = verbose)
   
   object <- spark.test(object = object,
                        check_positive = TRUE,
                        verbose = verbose)
   object@res_mtest
})


#' @param minimum.k minimum nearest neigbhours if maximum_distance != NULL, 
#' default, 2.
#' @param maximum.distance.delaunay distance cutoff for nearest neighbors 
#' to consider for Delaunay network, default 400
#' @param bin.method method to binarize gene expression, choice 
#' with kmeans and rank, defualt kmeans
#' @param expression.value expression values to use, choice with 
#' normalized, scaled and custom, defualt is normalized.
#' @param num.cores The number of cores to use, default 5.
#' @param python_Path python path. Path in Giotto object will be used when no extra dir is provided.
#' @param verbose verbosity of the function
#' @return A list contain theree methods of differentin expression
#' @importFrom Giotto createSpatialNetwork binSpect silhouetteRank
#' @importFrom reticulate import
#' @rdname DiffGenes
#' @method DiffGenes giotto

setMethod("DiffGenes", signature = "giotto",
          definition = function(object,
                                minimum.k = 2,
                                maximum.distance.delaunay = 400,
                                bin.method = "kmeans",
                                expression.value = "normalized",
                                num.cores = 5,
                                verbose = TRUE,
                                python_Path=NULL,
                                ...
                                ) {
   object <- createSpatialNetwork(gobject = object,
                                  minimum_k = minimum.k,
                                  maximum_distance_delaunay = maximum.distance.delaunay,
                                  ...)
   
   bin.method <- match.arg(bin.method, choices = c("kmeans", "rank"))
   expression.value <- match.arg(expression.value, choices = c("normalized", "scaled", "custom"))
   
   print("binspectGenes")
   binspectGenes <- binSpect(gobject = object, 
                              bin_method = bin.method,
                              expression_values = expression.value,
                              ...)
   
   print("silhouetteGenes")
   silhouetteGenes <- silhouetteRank(gobject = object,
                                     expression_values = expression.value,
                                     ...)
   print("spark")
   sparkObject <- GiottoToSpark(object = object)
   sparkGenes <- DiffGenes(object = sparkObject, 
                           num.cores = num.cores, 
                           ...)
   message("SpatialDE")
   count <- as.matrix(object@raw_exprs)
   count <- as.data.frame(count)
   locat <- as.data.frame(object@spatial_locs[, c("sdimx", "sdimy", "total_counts")])
   
   if(!is.null(python_Path)){
     reticulate::use_python(python = python_Path)
     py_config()}
   
   pd <- import("pandas")
   NaiveDE <- import("NaiveDE")
   SpatialDE <- import("SpatialDE")
   count <- count[rowSums(count) >= 20, ]
   norm_expr <- NaiveDE$stabilize(count)
   resid_expr <- NaiveDE$regress_out(locat, norm_expr, 'np.log(total_counts)')
   results <-  SpatialDE$run(locat, as.data.frame(t(resid_expr)))
   results <- results[order(results$qval), c("g", "l", "qval")]
   colnames(results)[1]="genes"
   
   DiffGene <- list(binspectGenes = binspectGenes,
                    silhouetteGenes = silhouetteGenes,
                    sparkGenes = sparkGenes,
                    spatialDE = results)
   DiffGene
})

#' @title gene pair coexpression based on scHOT
#' @param object a giotto object
#' @param fdrCutoff FDR value cutoff 
#' @param numGens number of genes for gene paire 
#' @param numCores The number of cores to use
#' @param ... other parameters for scHOT function \code{\link[scHOT]{scHOT}}
#' 
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom scater logNormCounts
#' @importFrom scran fitTrendVar modelGeneVar
#' @importFrom utils combn
#' @importFrom scHOT scHOT_buildFromSCE scHOT_addTestingScaffold scHOT_setWeightMatrix
#' @importFrom scHOT scHOT_calculateGlobalHigherOrderFunction scHOT_setPermutationScaffold
#' @importFrom scHOT scHOT_calculateHigherOrderTestStatistics scHOT_performPermutationTest
#' @importFrom scHOT scHOT
#' @importFrom BiocParallel MulticoreParam
#' 
#' @rdname CCI
#' @method CCI giotto

setMethod("CCI", signature = "giotto",
          definition = function(object,
                                fdrCutoff = 0.1,
                                numGens = 20,
                                numCores = 8,
                                ...) {
   position <- object@spatial_locs[, 1:2]
   colnames(position) <- c("x", "y")
             
   sce <- SingleCellExperiment(list(counts = as.matrix(object@raw_exprs)),colData = position)
   sce <- logNormCounts(sce)
   means <- rowMeans(logcounts(sce))
   vars <- rowVars(logcounts(sce))
   fit <- fitTrendVar(means, vars)
   var.out <- modelGeneVar(sce)
             
   seqvals <- seq(min(means), max(means), length.out = 1000)
   peakExp <- seqvals[which.max(fit$trend(seqvals))]
             
   hvg.out <- subset(var.out, subset = FDR <= fdrCutoff & mean > peakExp)  
   hvg.out <- hvg.out[order(hvg.out$bio, decreasing = TRUE), ]
             
   HVG <- sort(rownames(hvg.out))
   paris <- t(combn(HVG, 2))
   rownames(paris) <- paste(paris[, 1], paris[, 2], sep = "_")
             
   if (length(HVG) > numGens) {
      paris <- paris[sample(x = nrow(paris), size = 20), ]
   }
             
   scHOT_spatial <- scHOT_buildFromSCE(sce, assayName = "counts",
                                       positionType = "spatial",
                                       positionColData = c("x", "y"))
             
   scHOT_spatial <- scHOT_addTestingScaffold(scHOT_spatial, paris)
             
   scHOT_spatial <- scHOT_setWeightMatrix(scHOT_spatial,
                                          positionColData = c("x","y"),
                                          positionType = "spatial",
                                          nrow.out = NULL,
                                          span = 0.05)
   
   weightedPearson <- function (x, y, w = 1) {
      if (length(x) != length(y)) 
         stop("data must be the same length")
      if (length(w) == 1) {
         w <- rep(w, length(x))
      }
      nw = sum(w)
      wssx = nw * sum(w * (x^2)) - sum(w * x)^2
      wssy = nw * sum(w * (y^2)) - sum(w * y)^2
      wssxy = nw * sum(w * x * y) - sum(w * x) * sum(w * y)
      wcor = wssxy/sqrt(wssx * wssy)
      wcor
   }
   
   weightedSpearman <- function(x, y, w = 1) {
      if (length(x) != length(y)) {
         stop("x and y should have the same length")
      }
      if (length(w) == 1) {
         w <- rep(w, length(x))
      }
      keep = w > 0
      xr = rank(x[keep])
      yr = rank(y[keep])
      wp <- weightedPearson(x = xr, y = yr, w = w[keep])
      wp
   }
   
   
   
   scHOT_spatial <- scHOT_calculateGlobalHigherOrderFunction(scHOT_spatial, 
                                                             higherOrderFunction = weightedSpearman,
                                                             higherOrderFunctionType = "weighted")
             
   scHOT_spatial <- scHOT_setPermutationScaffold(scHOT_spatial, 
                                                 numberPermutations = 50,
                                                 numberScaffold = 10)
   scHOT_spatial <- scHOT_calculateHigherOrderTestStatistics(scHOT_spatial, 
                                                             higherOrderSummaryFunction = sd)
             
   BPPARAM = BiocParallel::MulticoreParam(workers = numCores)
   scHOT_spatial <- scHOT_performPermutationTest(scHOT_spatial, verbose = TRUE, 
                                                 parallel = TRUE, 
                                                 BPPARAM = BPPARAM)
             
   scHOT_spatial <- scHOT(scHOT_spatial, 
                          testingScaffold = paris,
                          positionType = "spatial",
                          positionColData = c("x", "y"),
                          nrow.out = NULL,
                          higherOrderFunction = weightedSpearman,
                          higherOrderFunctionType = "weighted",
                          numberPermutations = 50,
                          numberScaffold = 10,
                          higherOrderSummaryFunction = sd,
                          parallel = TRUE,
                          verbose = FALSE,
                          span = 0.05,
                          BPPARAM = BPPARAM,
                          ...)
             
   scHOT_spatial
})

#' @title Spatial domain visualization for cluster result.
#' @param object a giotto object
#' @param savePlot Boolean, whether to save plot.
#' @param outputFolder Output folder to save results. Path in Giotto object will be used when no external dir is provided.
#' @rdname DomainPlot
#' @import patchwork
#' @method DomainPlot giotto

setMethod("DomainPlot", signature = "giotto",
          definition = function(object,
                                savePlot=TRUE,
                                outputFolder=NULL,
                                ...
          ) {
            
            cell_metadata=as.data.frame(object@cell_metadata)
            g1=spotPOPplot(object=object,outputFolder=outputFolder, meta=cell_metadata$leiden,legend="leiden",title="Cell type annotation")
            g2=spotPOPplot(object=object,outputFolder=outputFolder, meta=cell_metadata$kmeans,legend="Kmeans",title="Spatial Clustering")
            hmrfindex=grep(pattern = "^hmrf",colnames(object@cell_metadata))[1]
            g3=spotPOPplot(object=object,outputFolder=outputFolder, meta=cell_metadata[,hmrfindex],legend="HMRF",title="Spatial Clustering")
            layout<-c(
              area(1,1,2,2),
              area(1,3,2,4),
              area(3,2,4,3))
            g4=( g2 + g3 + g1 )+ plot_layout(design=layout)
            g4=g4* theme_bw()
            if(savePlot){
              if(is.null(outputFolder)&!is.null(object@instructions$save_dir)){
                outputFolder=object@instructions$save_dir
              }
              
              ggsave(paste(outputFolder,"/DomainPlot.png",sep = ''),plot=g4,width=10,height=8,device=png,dpi=300)
            }
            return(g4)
          })






