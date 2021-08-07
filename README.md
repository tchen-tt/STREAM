

# **STREAM**: **S**patial **t**ranscriptomics data **a**nalysis and **m**odeling

**STREAM** is a protocol for processing spatial transcriptomic data. Currently, it supports
[10x-genomics](https://www.10xgenomics.com/products/spatial-gene-expression). The flow could also take processed matrices or Seurat object from
[FISH-seq](https://www.nature.com/articles/s41586-019-1049-y#Sec1), [STARmap](https://science.sciencemag.org/content/361/6400/eaat5691), [Slide-seq](https://science.sciencemag.org/content/363/6434/1463) and [Dbit-seq](https://www.cell.com/cell/pdf/S0092-8674(20)31390-8.pdf) [seq-scope](https://www.sciencedirect.com/science/article/pii/S0092867421006279?via%3Dihub) as input files. Processed matrices consist of count matrix and table with location information for each spots obtained from image processing. Besides，the meta file of cell-type annotations is optional.

SPADE is composed of for following parts: **Initialization**, **Preprocessing**, **SE genes**, **Spatial Clustering**,  **CCI** and **Integration with scRNA seq data**. `Preprocessing` part includes steps of initialization, visualization, QC, normalization, filter, dimension reduction, clustering based on expression only, marker-gene identification, etc. `SE genes` part is for the genes with strong difference in spatial expression pattern, in which published methods including `SpatialDE`, `binspect`, `silhouetteRank`, `SPARK` are collected. `Spatial clustering` part investigates sub-domains for tissue heterogeneity through `HMRF` and `Neighborhood-Kmeans method`. `CCI` part aims to study the potential gene pair relation and cell-cell interaction in the space.

## Installation


```{bash}
$ git clone https://github.com/YeehanXiao/STREAM.git
$ cd STREAM
$ conda env create -f environment.yml
$ conda activate STREAM
$ pip3 install NaiveDE
$ pip3 install spatialde
$ R -e 'devtools::install_github("xzhoulab/SPARK")'
$ R -e 'devtools::install_github("RubD/Giotto", force = TRUE)'
$ R -e 'devtools::install(".")'
```

## Usage

Initialize

``` {r1}
instrs = createGiottoInstructions(show_plot = FALSE,save_plot = FALSE, save_dir = "./",python_path = "/usr/local/bin/python3")

Giotto_obj <- createGiottoObject(raw_exprs = data[["Spatial"]]@counts, spatial_locs = data@meta.data[, c("x", "y")], instructions = instrs)
```

Preprocess

```{r2}
object <- processing(object = Giotto_obj)
```

Identify SE genes

```{r3}
diffgene <- DiffGenes(object = object, num.cores = 10)
```
Explore spatial clustering pattern

```{r4}
object <- KmeansCluster(object = object)
object <- HmrfCluster(object = object, outputFolder = "/Users/xiaoyihan", spatialGenes = diffgene$binspectGenes[1:200, ]$genes)
```
Study cell-cell interaction

```{r5}
cci <- CCI(object = object, numCores = 10)
interaction=SpatCellphoneDB(testob,num.cores = 16,cellphonedbPath = "/Users/xiaoyihan/miniconda3/bin/cellphonedb")
```

Visualization

```{r6}
QCplot(object=object)
SEplot(diffgene,Species = "mouse")
DomainPlot(object=object)
scHOTplot(cci)
cellphoneDBplot(interaction = interaction, outputFolder = "liver_results",save.plot = TRUE)
```

## References

  - [Dries, R., Zhu, Q. et al. Giotto: a toolbox for integrative analysis and visualization of spatial expression data. Genome Biol. 2021 Mar 8;22(1):78.](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02286-2)
  
  - [Zhu, Q., Dries, R. et al. Identification of spatially associated subpopulations by combining scRNAseq and sequential fluorescence in situ hybridization data. nature biotechnology. 2018 Oct 29.](https://www.nature.com/articles/nbt.4260)
  
  - [Valentine et al. SpatialDE: identification of spatially variable genes. Nature Methods. 2019 Mar 19.](https://www.nature.com/articles/nmeth.4636)
  
  ...
  
  
  
