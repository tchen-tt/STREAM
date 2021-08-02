#' Feature genens 
#' 
#' @param object An giotto object
#' @importFrom Giotto fDataDT
#' @export
Features <- function(object) {
  geneMetadata <- fDataDT(object)
  feature.gene <- subset(geneMetadata, hvg == "yes", select = gene_ID)
  feature.gene$gene_ID
}

#' Transform Giotto object to Spark
#'
#' @param object An Giotto object
#' @param ... Other arguments sell also \code{\link[SPARK]{CreateSPARKObject}}
#' 
#' @return An SPARK object
#' @importFrom  SPARK CreateSPARKObject
#'
#' @export

GiottoToSpark <- function(object, ...) {
  counts <- slot(object = object, name = "raw_exprs")
  position <- slot(object = object, name = "spatial_locs")[, 1:2]
  position <- as.data.frame(position)
  rownames(position) <- colnames(counts)
  spark <- CreateSPARKObject(counts = counts,
                             location = position,
                             ...)
  spark
}

#' @title TERM2GENE
#' @param category MSigDB collection abbreviation, such as H or C1.
#' @param Geneformat Format of gene name in count matrix, one of "ensembl_gene", "gene_symbol".
#' @param Species Species name, such as Homo sapiens or Mus musculus.
#' @importFrom   msigdbr msigdbr
#' @export term2gene
term2gene <- function(Geneformat=Geneformat,
                      Species=Species,
                      category=category){
  gene_sets = msigdbr(species = Species,category = category)
  if(Geneformat=="gene_symbol"){
    msigdbr_t2g = gene_sets %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()}
  if(Geneformat=="ensembl_gene"){
    msigdbr_t2g = gene_sets %>% dplyr::distinct(gs_name, entrez_gene) %>% as.data.frame()}
  return(msigdbr_t2g)
}




