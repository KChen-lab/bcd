library(R6)
library(Rtsne)
library(umap)
library(ggplot2)


# The core of varactor is distance.
# So there should be three fields:

rowSds <- function(x) sqrt(rowMeans((x - rowMeans(x)) ** 2) / dim(x)[2])

euclidean_pdist2 <- function(x, y){
  x <- as.matrix(x)
  if (missing(y)) y <- x
  else y <- as.matrix(y)
  xx = rowSums(x * x)
  yy = rowSums(y * y)
  pdist <- t(matrix(rep(yy, length(xx)), ncol=length(xx)))
  pdist <- pdist + xx
  pdist <- pdist - 2 * x %*% t(y)
  pdist
}

mahalanobis_pdist <- function(x, B){
  invB <- solve(B)
  xx = rowSums(x %*% invB * x)
  pdist <- t(matrix(rep(xx, length(xx)), ncol=length(xx)))
  pdist <- pdist + xx
  pdist <- pdist - 2 * x %*% invB %*% t(x)
  sqrt(pdist - min(pdist, 0))
}

categorical_boc_pdist <- function(expr, bad_label)
{
  expr <- as.matrix(expr)
  B <- matrix(0, ncol=dim(expr)[2], nrow=dim(expr)[2])
  for (j in unique(bad_label))
  {
    mask <- bad_label != j
    print(dim(expr))
    print(length(mask))
    temp <- t(expr[mask, ]) - colMeans(expr[!mask, ])
    B <- B + temp %*% t(temp)
  }
  B <- B / dim(expr)[1]
  mahalanobis_pdist(expr, B)
}

longitudinal_boc_pdist <- function(expr, bad_label, time_label, reg = 1., pow = 1.){
  B <- matrix(0, ncol=dim(expr)[2], nrow=dim(expr)[2])
  n <- dim(private$.reduced)[1]
  
  unique_bad_label <- unique(bad_label)
  time_of_bad_label <- vector(length=length(unique_bad_label))
  names(time_of_bad_label) <- unique_bad_label
  
  for (a_bad_label in unique_bad_label){
    correspond_time <- unique(time_label(bad_label == a_bad_label))
    if (length(correspond_time) == 1){
      time_of_bad_label[a_bad_label] <- correspond_time
    } else {
      stop('Times for cells in a batch must be the same!')
    }
  }
  
  W <- 1 / sqrt(euclidean_pdist2(x = time_of_bad_label, y = time_label) + reg) ** pow
  
  for (j in 1:length(unique_bad_label))
  {
    mask <- bad_label != unique_bad_label[j]
    temp <- t(expr[mask, ]) - colMeans(expr[!mask, ])
    temp <- t(t(temp) * W[j, mask])
    B <- B + temp %*% t(temp)
  }
  B <- B / sum(W ** 2)
  mahalanobis_pdist(expr, B)
}

RunBoc <- function(object, batch, method='categorical', from='pca', to='umap', 
                   reduction.name='boc', day, verbose=TRUE){
  if (verbose) cat("Calculating pairwise distances...")
  bad_labels <- lapply(X = object[[batch]], function(x) sub(x, ' ', '_'))
  bad_label <- do.call(paste, unname(object[[batch]]))
  if (method == 'categorical'){
    pdist <- categorical_boc_pdist(object@reductions[[from]]@cell.embeddings, bad_label)
  } else if (method == 'longitudinal') {
    if(missing(day)){
      stop('For longitudinal data, parameter "day" is required')
    }
    pdist <- longitudinal_boc_pdist(object@reductions[[from]], bad_label, object[[day]])
  } else {
    stop('Unknown distance. Please choose from categorical and longitudinal.')
  }
  
  if (verbose) cat("Calculating embedding...")
  object@reductions[['boc']] <- new('DimReduc')
  object@reductions[['boc']]@key <- 'BOC_'
  object@reductions[['boc']]@assay.used <- object@reductions[[from]]@assay.used
  if (to == 'umap'){
    config <- umap.defaults
    config$input <- "dist"
    object@reductions[['boc']]@cell.embeddings = umap(pdist, config = config)$layout
  } else if (to == 'umap3d'){
      config <- umap.defaults
      config$input <- "dist"
      config$n_components <- 3
      object@reductions[['boc']]@cell.embeddings = umap(pdist, config = config)$layout
  } else if (to == 'tsne'){
    object@reductions[['boc']]@cell.embeddings = Rtsne(pdist, is_distance = T)$Y
  } else {
    stop('Unknown embedding method. Please choose from umap, umap3d and tsne.')
  }
  return(object)
}
