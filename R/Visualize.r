#' @title spotPOPplot
#' @param object A giotto object
#' @param outputFolder Output folder to save results
#' @param title A character for the title of plot.
#' @param meta A discrete vector to assign color for each spot, e.g., cell type annotation, HMRF result, expression of a certain gene.
#' @param legend A character of feature to divide up cell into sub-populations, e.g. "K-means", "HMRF", "cell type".
#' @param savePlot Boolean, whether to save plot.
#' @param ... other arguments 
#' @import patchwork
#' @return ggplot object
#' @export
#' @rdname Visualize

spotPOPplot <- function(object, 
                        outputFolder = NULL, 
                        title = "CellPOP",
                        meta = object@cell_metadata$leiden, 
                        legend = "leiden", 
                        savePlot = FALSE, 
                        ...){
   Nspots <- ncol(object@raw_exprs)
   g <- ggplot()
   g <- g + geom_point(aes(x=object@spatial_locs$sdimx,y=object@spatial_locs$sdimy,color=as.factor(meta)),pch = 20, size = 1500/Nspots)
   g <- g + labs(x = "x", y = "y", title=title, color=legend)
   g <- g + scale_color_npg(alpha=0.7)
   g <- g * theme_bw()
   if(savePlot){
      ggsave(paste0(outputFolder,"/",title,".png"),plot = g, device = NULL,width = 7)
      }
   return(g)
}

#' @title Scoreplot
#' @description Visualize the spatial data with a continuous value in each spot, e.g. expression, scHOT score.
#' @param x Numeric, x-coordinate of spatial data.
#' @param y Numeric, y-coordinate of spatial data.
#' @param outputFolder Output folder to save results
#' @param title A character for the title of plot.
#' @param value A vector to assign color for each spot, e.g., cell type annotation, HMRF result, expression of a certain gene.
#' @param legend A character of feature to divide up cell into sub-populations, e.g. "K-means", "HMRF", "cell type".
#' @param savePlot Boolean, whether to save plot.
#' @param ... other arguments 
#' @import ggsci
#' @return ggplot object
#' @export
#' @rdname Visualize

Scoreplot=function(x,y,outputFolder=NULL,title="Spatplot",value,legend="score",savePlot=FALSE,...){
    Nspots=length(x) 
    p<-ggplot()
    p <- p + geom_point(aes(x=x,y=y,color=value), pch = 20, size = 1500/Nspots)
    p <- p + labs(x = "x", y = "y", title=title)
    p <- p + scale_color_gradientn(name=legend,colours = c(pal_simpsons("springfield")(10)[c(3,2,10,6,1,5)]))
    p <- p * theme_bw()
    if(savePlot){
       ggsave(paste0(outputFolder,"/",title,".png"),plot = g, device = NULL,width = 7)
    }
    return(p)
 }
 
#' @title QCplot
#' @param object A giotto object
#' @param outputFolder Output folder to save results
#' @param ... other arguments
#' @export
#' @rdname Visualize
 
QCplot=function(object,outputFolder=NULL,...){
   if(is.null(outputFolder)&!is.null(object@instructions$save_dir)){
      outputFolder=object@instructions$save_dir
   }
    spotsvalue=log(colSums(object@raw_exprs),base = 10)
    genesvalue=log(apply(object@raw_exprs,2,function(x){length(which(x!=0))}),base =10)
    mitoindex=grep("mt-",rownames(object@raw_exprs))
    if(length(k)<5) {
      print("Little mito genes detected")
      } else{
       mitopercent=colSums(object@raw_exprs[mitoindex,])/colSums(object@raw_exprs)
       g <- ggplot()
       g <- g + geom_point(aes(x=object@spatial_locs$sdimx,y=object@spatial_locs$sdimy,color=mitopercent),pch =20,size = 1500/Nspots)
       g <- g + labs(x = "x", y = "y", title="MT reads fraction")
       g <- g + scale_color_gradientn(name = "fraction",colours = c(pal_simpsons("springfield")(10)[c(3,2,10,6,1,5)]))
       g <- g * theme_bw()
       ggsave(paste0(outputFolder,"/MTfraction.png"),plot = g, device = NULL,width = 5,dpi=300)
    }
    Nspots=ncol(object@raw_exprs)
    g1 <- ggplot()
    g1 <- g1+ geom_point(aes(x=object@spatial_locs$sdimx,y=object@spatial_locs$sdimy,color=spotsvalue),pch = 20, size = 1500/Nspots)
    g1 <- g1 + labs(x = "x", y = "y", title="Counts per spot")
    g1 <- g1 + scale_color_gradientn(name="log 10 counts\n ",colours = c(pal_simpsons("springfield")(10)[c(3,2,10,6,1,5)]))
    
    g2 <- ggplot()
    g2 <- g2 + geom_point(aes(x=object@spatial_locs$sdimx,y=object@spatial_locs$sdimy,color=genesvalue),pch = 20,size = 1500/Nspots)
    g2 <- g2 + labs(x = "x", y = "y", title="Genes detected per spots")
    g2 <- g2 + scale_color_gradientn(name = "log 10 genes",colours = c(pal_simpsons("springfield")(10)[c(3,2,10,6,1,5)]))
    
    
    g3 <- (g1 + g2) * theme_bw()
    ggsave(paste0(outputFolder,"/QC_plot.png"),plot = g3, device = NULL,width = 10,height=4,dpi=300)
 }
 
#' @title scHOTplot
#' @param scHOTobj A scHOT object
#' @param outputFolder Output folder to save results. 
#' @param genes A gene pair in scHOT object with weighted higher order statistc, linked by "_" or taken as a vector.
#' @param savePlot Boolean, whether to save plot.
#' @param ... other arguments 
#' @import  ggsci
#' @importFrom  cowplot plot_grid
#' @return ggplot object
#' @export
#' @rdname Visualize
scHOTplot=function(scHOTobj,outputFolder=NULL,title="gene-gene coexpression",savePlot=FALSE,genes=NULL,...){
    out<- scHOTobj@scHOT_output[order(scHOTobj@scHOT_output$pvalEstimated, decreasing=TRUE),] 
    if(is.null(genes)){
       genes=rownames(out)[1]}
    if (!("higherOrderSequence" %in% colnames(scHOTobj@scHOT_output))) {
       stop("higherOrderSequence is not found in scHOT_output")
    }
    namescols = grep("gene_", colnames(scHOTobj@scHOT_output), value = TRUE)
    wcor <- as.matrix(scHOTobj@scHOT_output$higherOrderSequence)
    rownames(wcor) <- apply(scHOTobj@scHOT_output[, namescols, drop = FALSE], 
                            1, paste0, collapse = "_")
    colnames(wcor) <- NULL
    if (ncol(wcor) != ncol(scHOTobj)) {
       warning("Not all the cell position has higherOrderSequence statistics,\n            set nrow.out = NULL in scHOT_setWeightMatrix to calculate\n            higherOrderSequence for all positions!")
    }
    else {
       if (!genes[1] %in% rownames(wcor)) {
          if (length(genes) == 2) {
             if (!paste0(sort(genes), collapse = "_") %in% 
                 rownames(wcor)) {
                stop("gene pair has no higherOrderSequence ")
             }else {
                gene1=genes[1]
                gene2=genes[2]
                genes <- paste0(sort(genes), collapse = "_")
             }
          }else {
             stop("gene pairs has no higherOrderSequence ")
          }
       }else{
          gene1=scHOTobj@scHOT_output[genes, namescols[1]]
          gene2=scHOTobj@scHOT_output[genes, namescols[2]]
       }}
    wcor_all <- wcor[genes, , drop = FALSE]
    p1=Scoreplot(x = scHOTobj@colData$x ,y = scHOTobj@colData$y,outputFolder = outputFolder,title = genes,value = wcor_all,legend = "scHOT score",savePlot =FALSE)
    #gene1
    vgene1 <- scHOTobj@assays@data[[1]][gene1,]
    vgene1.q95 <- quantile(vgene1, 0.95)
    vgene1[vgene1 > vgene1.q95] <- vgene1.q95
    p2=Scoreplot(x = scHOTobj@colData$x ,y = scHOTobj@colData$y,outputFolder = outputFolder,value=vgene1,legend="Expression",title=gene1,savePlot = FALSE)
    #gene2
    vgene2 <- scHOTobj@assays@data[[1]][gene2,]
    vgene2.q95 <- quantile(vgene2, 0.95)
    vgene2[vgene2 > vgene2.q95] <- vgene2.q95
    p3=Scoreplot(x = scHOTobj@colData$x ,y = scHOTobj@colData$y,outputFolder = outputFolder,value=vgene2,legend="Expression",title=gene2,savePlot = FALSE)
    p = plot_grid(p1, NULL, p2,p3, ncol = 2)
    p=p*theme_bw() * theme(panel.border = element_blank())
    if(savePlot){
       ggsave(paste0(outputFolder,"/scHOT_",genes,".png"),plot = p, device = png,width = 10)
    }
    return(p)
}



#' @title SEgenesplot
#' @param diffgenes A list out of SE-genes calculation.
#' @param topgenes Integer; set how many top SE-genes to take, default as 100.
#' @param method A character to select the SE method for downstream analysis, one of "Binspect","silhouetteRank","SPARK","Integrated".
#' @param outputFolder Output folder to save results.
#' @param savePlot Boolean, whether to save plot.
#' @param category MSigDB collection abbreviation, such as H or C1.
#' @param Geneformat Format of gene name in count matrix, one of "ensembl_gene", "gene_symbol".
#' @param Species Species name, such as Homo sapiens or Mus musculus.
#' @param ... other arguments 
#' @import ggVennDiagram
#' @importFrom clusterProfiler enricher 
#' @import patchwork
#' @import parallel
#' @return ggplot object
#' @export SEplot
#' @rdname Visualize
SEplot=function(diffgenes=diffgenes, 
                 topgenes=100,
                 method="Integrated", 
                 outputFolder=NULL, 
                 savePlot=TRUE,
                 Species="human",
                 category="C2", 
                 Geneformat="gene_symbol",
                 ...){
  
  if(is.null(outputFolder) & !is.null(object@instructions$save_dir)){
    outputFolder=object@instructions$save_dir}
    
   method <- match.arg(method, choices = c("Binspect","silhouetteRank","SPARK","Integrated")) 
   if(method=="Integrated"){
     sparkres <- diffgenes[["sparkGenes"]]
     spark_spatialgenes<-rownames(sparkres[order(sparkres$adjusted_pvalue,decreasing = FALSE),])
     topspatialDEgenes<-diffgenes[["spatialDE"]]$genes[1:topgenes]#1
     topbingenes<-diffgenes[["binspectGenes"]]$genes[1:topgenes] #3
     topsilgenes<-diffgenes[["silhouetteGenes"]]$genes[1:topgenes]#2
     topsparkgenes<-spark_spatialgenes[1:topgenes]#4
     X=list(
       spatialDE=topspatialDEgenes,
       silhouetteRank=topsilgenes,
       Binspect=topbingenes,
       Spark=topsparkgenes)
     vn = ggVennDiagram(X, label_alpha = 0,label_color = "white")
     vn = vn + theme_bw() + theme(panel.border = element_blank(),panel.grid = element_blank(), 
                                  axis.text = element_blank(), 
                                  axis.ticks = element_blank())
     
     #filter spatialgenes
     ALL<-Reduce(union,X)
     sum_matrix <- matrix(data=0,nrow=length(ALL),ncol=4)
     rownames(sum_matrix)<-ALL
     index=lapply(c(seq(1,4)),function(x){
       sum_matrix[as.vector(unlist(X[x])),x]=1
       return(sum_matrix[,x])
     })
     sum_matrix<-do.call(cbind,index)
     colnames(sum_matrix)<- c("spatialDE","silRank", "binspect","spark")
     sums=rowSums(sum_matrix)
     core_spatialgenes<-sort(sums,decreasing = TRUE)[1:4]
     final_spatialgenes<-names(sort(sums,decreasing = TRUE)[1:topgenes])
     
     spatlist=mclapply(names(core_spatialgenes),function(gene){
       p=Scoreplot(x = object@spatial_locs$sdimx ,y = object@spatial_locs$sdimy,outputFolder = outputFolder,title = gene,value = log2(object@raw_exprs[gene,]+1),legend = "log2(Expression)",savePlot =FALSE)
       return(p)})
     if(savePlot){
       ggsave(paste(outputFolder,"/vnPlot.png",sep = ''),plot=vn,width=10,device=png,dpi=300,scale=0.8)
     }
     
   }
   if(method=="Binspect"){
     if(!length(diffgenes[["binspectGenes"]])){stop("No Binspect list in diffgenes")}
     spatlist=mclapply(diffgenes[["binspectGenes"]]$gene[1:4],function(gene){
       p=Scoreplot(x = object@spatial_locs$sdimx ,y = object@spatial_locs$sdimy,outputFolder = outputFolder,title = gene,value = log2(object@raw_exprs[gene,]+1),legend = "log2(Expression)",savePlot =FALSE)
       return(p)})
   }
   
   if(method=="silhouetteRank"){
     if(!length(diffgenes[["silhouetteGenes"]])){stop("No silhouette list in diffgenes")}
     spatlist=mclapply(diffgenes[["silhouetteGenes"]]$genes[1:4],function(gene){
       p=Scoreplot(x = object@spatial_locs$sdimx ,y = object@spatial_locs$sdimy,outputFolder = outputFolder,title = gene,value = log2(object@raw_exprs[gene,]+1),legend = "log2(Expression)",savePlot =FALSE)
       return(p)})
   }
   
   if(method=="SPARK"){
     if(!length(diffgenes[["sparkGenes"]])){stop("No silhouette list in diffgenes")}
     sparkres <- diffgenes[["sparkGenes"]]
     spark_spatialgenes<-rownames(sparkres[order(sparkres$adjusted_pvalue,decreasing = FALSE),])
     spatlist=mclapply(spark_spatialgenes$genes[1:4], function(gene){
       p=Scoreplot(x = object@spatial_locs$sdimx ,y = object@spatial_locs$sdimy,outputFolder = outputFolder,title = gene,value = log2(object@raw_exprs[gene,]+1),legend = "log2(Expression)",savePlot =FALSE)
       return(p)})
   }
   
   if(method=="spatialDE"){
     if(!length(diffgenes[["spatialDE"]])){stop("No silhouette list in diffgenes")}
     spatlist=mclapply(diffgenes[["spatialDE"]]$genes[1:4], function(gene){
       p=Scoreplot(x = object@spatial_locs$sdimx ,y = object@spatial_locs$sdimy,outputFolder = outputFolder,title = gene,value = log2(object@raw_exprs[gene,]+1),legend = "log2(Expression)",savePlot =FALSE)
       return(p)})
   }
   
   gene_sets = msigdbr(species = Species,category = category)
   
   if(Geneformat=="gene_symbol"){
     msigdbr_t2g = gene_sets %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()}
   if(Geneformat=="ensembl_gene"){
     msigdbr_t2g = gene_sets %>% dplyr::distinct(gs_name, entrez_gene) %>% as.data.frame()}
   
   enrich <- enricher(gene = final_spatialgenes, TERM2GENE = msigdbr_t2g)
   
   SE <- plot_grid(spatlist[[1]],spatlist[[2]],spatlist[[3]],spatlist[[4]],ncol=2)
   SEenrichdot <- clusterProfiler::dotplot(enrich,x="GeneRatio",color="p.adjust",showCategory=10,
                       size=NULL,split=NULL,font.size=12,title="KEGG of SEgenes")
   SEenrichdot <- SEenrichdot + theme(axis.text.y = element_text(size=6)) + 
      scale_x_continuous(expand = c(0, 0.1))
   if(savePlot){
     ggsave(paste(outputFolder,"/enrichDot.png",sep = ''),plot=SEenrichdot)
     ggsave(paste(outputFolder,"/SEgenes_",method,".png",sep = ''),plot=SE,width=10,height=8,device=png,dpi=300)
   }
}

#' @title cellphoneDBplot
#' @param interaction A list out of cellphoneDB calculation.
#' @param Lcelltype 
#' @param Rcelltype 
#' @param top Number of top pairs to plot.
#' @param outputFolder Output folder to save results.
#' @param save.plot directly save the plot [boolean].
#' @param ... other arguments 
#' @importFrom  pheatmap pheatmap
#' @export cellphoneDBplot
#' @rdname Visualize
cellphoneDBplot = function(interaction=interaction, outputFolder=NULL, save.plot=FALSE,...){
  
  #heatmap
  cellNetwork <- read.delim(interaction$countNetwork, header = TRUE, stringsAsFactors = FALSE)
  cellNetwork.spread <-  cellNetwork %>% spread(key = TARGET, value = count)
  rownames(cellNetwork.spread) <- cellNetwork.spread$SOURCE
  cellNetwork.spread <- cellNetwork.spread[, !(colnames(cellNetwork.spread) %in% "SOURCE")]
  if (save.plot){
    png(filename=paste0(outputFolder,"/SpatcellphoneDBheatmap.png"))
    pheatmap::pheatmap(cellNetwork.spread, angle_col=0)
    dev.off()
  }else{pheatmap::pheatmap(cellNetwork.spread, angle_col=0)}
  #dotplot
  p2=cellphondedb_dotplot(Lcelltype = 1, Rcelltype = 1:4, interaction = interaction,top = 5)
  p2
  if (save.plot){
    ggsave(filename=paste0(outputFolder,"/SpatcellphoneDBdot.png"),plot = p2)
  }
}




#' @title basic cellphonedb dotplot
#' @param interaction A list out of cellphoneDB calculation.
#' @param selected_columns selected genepairs of cellphoneDB results
#' @param selected_rows selected clusters of cellphoneDB results
#' @param means_separator separator of means value file
#' @param pvalues_separator separator of p-value file
#' @param ... other arguments 
#' @import ggplot2
#' @importFrom grDevices colorRampPalette
#' @rdname Visualize

dot_plot = function(
                    interaction,
                    selected_rows = NULL,
                    selected_columns = NULL,
                    means_separator = '\t',
                    pvalues_separator = '\t',
                    ...
){
  
  means_path = interaction$means
  pvalues_path = interaction$pvalue
  all_pval = read.delim(pvalues_path, header=T, stringsAsFactors = F, sep=means_separator, comment.char = '', check.names=F)
  all_means = read.delim(means_path, header=T, stringsAsFactors = F, sep=pvalues_separator, comment.char = '', check.names=F)
  
  intr_pairs = all_pval$interacting_pair
  all_pval = all_pval[,-c(1:11)]
  all_means = all_means[,-c(1:11)]
  
  if(is.null(selected_rows)){
    selected_rows = intr_pairs
  }
  
  if(is.null(selected_columns)){
    selected_columns = colnames(all_pval)
  }
  
  sel_pval = all_pval[match(selected_rows, intr_pairs), selected_columns]
  sel_means = all_means[match(selected_rows, intr_pairs), selected_columns]
  
  df_names = expand.grid(selected_rows, selected_columns)#combination
  pval = unlist(sel_pval)
  pval[pval==0] = 0.0009
  plot.data = cbind(df_names,pval)
  pr = unlist(as.data.frame(sel_means))
  pr[pr==0] = 1
  plot.data = cbind(plot.data,log2(pr))
  colnames(plot.data) = c('pair', 'clusters', 'pvalue', 'mean')
  
  my_palette <- colorRampPalette(c("black", "blue", "yellow", "red"), alpha=TRUE)(n=399)
  
  dotplot=ggplot(plot.data,aes(x=clusters,y=pair)) +
    geom_point(aes(size=-log10(pvalue),color=mean)) +
    scale_color_gradientn('Log2 mean (Molecule1, Molecule2)', colors=my_palette) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text=element_text(size=14, colour = "black"),
          axis.text.x = element_text(angle = 90, hjust = 1),
          axis.text.y = element_text(size=12, colour = "black"),
          axis.title=element_blank(),
          panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))
  return(dotplot)
}


#' @title cellphondedb_dotplot
#' @param interaction A list out of cellphoneDB calculation.
#' @param Lcelltype ligand cell type
#' @param Rcelltype Receptor cell type
#' @param top Number of top pairs to plot.
#' @param outputFolder Output folder to save results.
#' @param save.plot directly save the plot [boolean]
#' @param ... other arguments 
#' @importFrom  pheatmap pheatmap
#' @rdname Visualize
cellphondedb_dotplot <- function(Lcelltype = 1, Rcelltype = 2, interaction, top = 10) {
  Rcelltype.Lcelltype <- paste(Rcelltype, Lcelltype, sep = "|")
  Lcelltype.Rcelltype <- paste(Lcelltype, Rcelltype, sep = "|")
  cellphonedb.merged<-cellphonedb_merge(interaction)
  cellphone <- na.omit(cellphonedb.merged)
  Rcelltype.Lcelltype.top <- cellphone[cellphone$Interaction %in% Rcelltype.Lcelltype, ]
  if (nrow(Rcelltype.Lcelltype.top) >= top)
    Rcelltype.Lcelltype.interaction <- Rcelltype.Lcelltype.top[order(Rcelltype.Lcelltype.top$mean, decreasing = TRUE)[1:top], "interacting_pair", drop = TRUE]
  else
    Rcelltype.Lcelltype.interaction <- Rcelltype.Lcelltype.top[, "interacting_pair", drop = TRUE]
  Lcelltype.Rcelltype.top <- cellphone[cellphone$Interaction %in% Lcelltype.Rcelltype, ]
  if (nrow(Lcelltype.Rcelltype.top) >= top)
    Lcelltype.Rcelltype.interaction <- Lcelltype.Rcelltype.top[order(Lcelltype.Rcelltype.top$mean, decreasing = TRUE)[1:top], "interacting_pair", drop = TRUE]
  else
    Lcelltype.Rcelltype.interaction <- Lcelltype.Rcelltype.top[, "interacting_pair", drop = TRUE]
  dot_plot(selected_rows = c(Rcelltype.Lcelltype.interaction, Lcelltype.Rcelltype.interaction), 
           selected_columns = c(Rcelltype.Lcelltype, Lcelltype.Rcelltype), 
           interaction=interaction)
}

