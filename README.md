# Introduction
Varactor (Visualization and analysis of single-cell RNA-seq data by alternative clustering) takes the expression matrix and the alleged nonpreferred clustering, which are either labeled beforehand, or previously discovered from the data. It redefines the pairwise distance of the cells, where the effect of the unwanted clustering is controlled, and use tSNE (which can be substituted with any distance-based clustering methods) on the new distances. It provides a new perspective for inspecting the similarity of the cells, instead of manipulating the expression data. The differential expression analysis may then be performed on the original data (see Nygaard, V. et al. Methods that remove batch effects while retaining group differences may lead to exaggerated confidence in downstream analyses. Biostatistics, 17(1), 29–39 (2016)) with the batches considered either covariates or strata. 

In the case one only wants to use the stratified test, we implemented stratified Wilcoxon U-test as a supplemnent to the [Seurat](https://github.com/satijalab/seurat) package. Please see the corresponding [repository](https://github.com/KChen-lab/stratified-tests-for-seurat).

# Usage
This is an abstract of the analysis of the PBMC datasets. The full examples is available [here](https://kchen-lab.github.io/varactor/pbmc_varactor.nb.html).

## Create a Varactor Object
Varactor is implemented using R6 class. To create an object, simply call ```new()```. The following chunk uses 10x and well-seq of PBMC for example, to create a Varactor object. A label called "sample" is automatically added to the labels to flag the origination of cells, even after they are combined into one matrix (by ```combine()``` showing later).
```r
source("./R/varactor_class.R")
data <- list(x10x = pbmc_10x$expr, well = pbmc_well$expr)
labels <- list(x10x = list(type = as.character(pbmc_10x$cell_type)), 
                           well = list(type = as.character(pbmc_well$cell_type))
                          )
obj2 <- Varactor$new(data = data, labels = labels)
```

## Preprocess
The input is a list of datasets and labels, so preprocessing is necessary to analyze them jointly. This includes normalization (including log-transform) of each dataset, combine them together (while only retain genes appear in all datasets) and dimensional reduction. This can be done by calling the corresponding methods (i.e., member functions).
```r
obj$normalize()
obj$combine()
obj$reduce()
```
Altenatively, you may chain these methods as follows.
```r
obj$normalize()$combine()$reduce()
```
The object ```obj``` will contain the processed dataset and you do not need to assign the value back as in ~~```obj <- obj$normalize()```~~.

## Embedding and clustering
All process are routine till now. To stress the idea of "alternative clustering", you can now create branches in the object to embed and cluster data with different metrics. 
### Primary
Using Euclidean distance, it is easy to find and plot the most familiar embedding using t-SNE. You may also choose to use UMAP by using ```embed("primary", "umap")```.
```r
obj2$define_metric("primary", "euclidean")$measure("primary")$embed("primary", "tsne")$plot_embedding("primary", "type", pch=20)
```
You can also plot it with a different coloring, to emphasize the difference between samples instead of cell types. The embedding is automatically stored in the object, so you only need to call ```plot_embedding```.
```{r}
obj2$plot_embedding("primary", "sample", pch=20)
```

### Alternative
By using a different definition of distance (davidson distance controling difference between samples in this case), you can find the alternative embedding.
```r
obj2$define_metric("alternative", "davidson", strata = "sample")$measure("alternative")$embed("alternative", "tsne")$plot_embedding("alternative", "type", pch=20)
```

## Hypothesis Testing for Differential Expression Genes
We provide two types of hypothesis tests, the ones based on strata and the ones based on covariates.

## Downstream Trajectory Inference


## Non-distance-based methods?
The metric based nature of Varactor makes it suitable for distance based methods. That said, there are workarounds for it to fit non-distance-based methods. The easiest one is to perform a multidimensional scaling (MDS) to transfer the distance matrix back in to a space.

# Full Examples

## PBMC cells
This [experiment](https://kchen-lab.github.io/varactor/pbmc_varactor.nb.html) show that Varactor can account for samples assayed by different technologies.

[10x PBMC3k dataset](http://support.10xgenomics.com/single-cell/datasets/pbmc3k)

[Seq-Well PBMC dataset (GSE92495)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92495)

## Brain cells
This experiment show that Varactor can evoke findings on multi-sample data.

Lake, B. B. et al. Neuronal subtypes and diversity revealed by single-nucleus RNA sequencing of the human brain. Science 352, 1586–1590 (2016)

Many other good datasets are available at [Hemberg Group, Sanger Institute](https://github.com/hemberg-lab/scRNA.seq.datasets).
