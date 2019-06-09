# varactor
Visualization and analysis of single-cell RNA-seq data by alternative clustering

# Introduction
Varactor takes the expression matrix and the alleged nonpreferred clustering, which are either labeled beforehand, or previously discovered from the data. It redefines the pairwise distance of the cells, where the effect of the unwanted clustering is controlled, and use tSNE (which can be substituted with any distance-based clustering methods) on the new distances. It provides a new perspective for inspecting the similarity of the cells, instead of manipulating the expression data. The differential expression analysis may then be performed on the original data (see Nygaard, V. et al. Methods that remove batch effects while retaining group differences may lead to exaggerated confidence in downstream analyses. Biostatistics, 17(1), 29–39 (2016)) with the batches considered either covariates or strata. We implemented stratified Wilcoxon U-test as a stand alone tool, and a supplemnent to the [Seurat](https://github.com/satijalab/seurat) package.

# Results
## Datasets
Many good datasets can be found from [Hemberg Group, Sanger Institute](https://github.com/hemberg-lab/scRNA.seq.datasets)

## Brain cells
Lake, B. B. et al. Neuronal subtypes and diversity revealed by single-nucleus RNA sequencing of the human brain. Science 352, 1586–1590 (2016)
