#' @title cellphoneDB for spot-population
#' @param object A giotto object
#' @param outputFolder Output folder to save results.
#' @param num.cores The number of cores to use, default as 4.
#' @param cellphonedbPath Path to launch cellphoneDB.
#' @param Geneformat Format of gene name in count matrix, one of "ensembl_gene", "gene_symbol".
#' @param method Spatial clustering method to annotate spots population, either "kmeans" or "HMRF".
#' @param ... other arguments 
#' @export
#' @rdname SpatCCI

SpatCellphoneDB <- function(object, 
                            outputFolder=NULL,
                            num.cores = 4,
                            Geneformat = "hgnc_symbol",
                            cellphonedbPath = NULL,
                            method = "kmeans",
                            ...) {
  method <- match.arg(method, c("kmeans","HMRF","leiden"))
  Geneformat <- match.arg(Geneformat, c("ensembl", "gene_name", "hgnc_symbol"))
  transformcounts <- object@raw_exprs
  transformcounts <- as.data.frame(as.matrix(transformcounts))
  colnames(transformcounts) <- paste0("c", colnames(transformcounts))
  colnames(transformcounts) <- gsub(pattern = "_", replacement = "", x = colnames(transformcounts))
  
  counts <- cbind.data.frame(data.frame(Gene = toupper(rownames(transformcounts))),
                           transformcounts)
  
  meta <-  data.frame(Cell=object@cell_metadata$cell_ID)
  cell_metadata <- as.data.frame(object@cell_metadata)
  if (method == "kmeans") {
    kmeans <- grep("kmeans", colnames(object@cell_metadata))[1]
    if (is.null(kmeans)) {
      stop("kmeans is not contains in cell_meatadata")
    }
    meta$cell_type <- object@cell_metadata[, kmeans]
  } else if (method=="HMRF"){
    hmrfindex <- grep(pattern = "^hmrf", colnames(object@cell_metadata))[1]
    if (is.null(hmrfindex)) {
      stop("hmrf is not cantains in cell_metadata")
  }
    meta$cell_type <- cell_metadata[, hmrfindex]
  } else if (method=="leiden"){
      leiden <- grep(pattern = "leiden", colnames(object@cell_metadata))[1]
      if (is.null(leiden)) {
        stop("leiden is not cantains in cell_metadata")
      }
    meta$cell_type <- cell_metadata[, leiden]
  }
  meta$Cell <- paste0("c", meta$Cell)
  meta$Cell <- gsub(pattern = "_", replacement = "", x = meta$Cell)
  
  counts <- counts[, c("Gene", meta$Cell)]
  
  if (is.null(outputFolder) & !is.null(object@instructions$save_dir)) {
    outputFolder=object@instructions$save_dir
  }
  outputFolder <- paste0(outputFolder,"/cellphonedb")
  
  if (!dir.exists(outputFolder)) {
    dir.create(outputFolder, recursive = TRUE)
  }
  
  write.table(counts, file = paste0(outputFolder,"/counts.txt"), sep= "\t", quote = FALSE, row.names = FALSE)
  write.table(meta, file = paste0(outputFolder,"/meta.txt"), sep= '\t', quote = FALSE, row.names = FALSE)
  
  commands <- paste0(cellphonedbPath, " ")
  methods <- "method statistical_analysis "
  gene <- paste0("--counts-data ", Geneformat, " ")
  
  count <- paste0(outputFolder, "/counts.txt ")
  metads <- paste0(outputFolder, "/meta.txt ")
  
  thread <- paste0("--threads ", num.cores, " ")
  output <- paste0("--output-path ", outputFolder)
  system(paste0(commands, methods, gene, metads, count, thread, output))
  system(paste0(commands, "plot heatmap_plot ", "--pvalues-path ", paste0(outputFolder, "/pvalues.txt", " "), 
                metads, " --output-path ", outputFolder))
  list(countNetwork = paste0(outputFolder, "/count_network.txt"),
       pvalue = paste0(outputFolder, "/pvalues.txt"),
       singnificantMeans = paste0(outputFolder, "/significant_means.txt"),
       means = paste0(outputFolder, "/means.txt"))
}

#' @title cellphonedb_merge
#' @param interaction A list out of cellphoneDB calculation.
#' @param ... other arguments 
#' @importFrom magrittr `%>%`
#' @rdname SpatCCI

cellphonedb_merge=function(interaction){
  cellphone.pvalues <- read.delim(interaction$pvalue, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  cellphone.pvalues.spread <- cellphone.pvalues %>%
    gather(-(id_cp_interaction:is_integrin), key = "Interaction", value = "pvalue") %>%
    select(interacting_pair, Interaction, pvalue)
  
  cellphone.sigMeans <- read.delim(interaction$singnificantMeans, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  cellphone.sigMeans.spread <- cellphone.sigMeans %>%
    gather(-(id_cp_interaction:rank), key = "Interaction", value = "mean") %>%
    select(interacting_pair, Interaction, mean)
  
  cellphone.pvalues.spread$mean <- cellphone.sigMeans.spread$mean[match(paste0(cellphone.pvalues.spread$interacting_pair, cellphone.pvalues.spread$Interaction), paste0(cellphone.sigMeans.spread$interacting_pair, cellphone.sigMeans.spread$Interaction))]
  return(cellphone.pvalues.spread)
}

