
RunALT <- function(object, dims.use, feature.unwanted, feature.time.dict, reduction.use, reduction.name = 'alt', 
                   reduction.key = 'ALT_', reg=1, pow=1, strength=1) {
  euclidean_pdist2 <- function(x, y){
    x <- as.matrix(x)
    if (missing(y)) {
      y <- x
    } else {
      y <- as.matrix(y)
    } 
    xx = rowSums(x * x)
    yy = rowSums(y * y)
    pdist2 <- t(matrix(rep(yy, length(xx)), ncol=length(xx)))
    pdist2 <- pdist2 + xx
    pdist2 <- pdist2 - 2 * x %*% t(y)
    pdist2
  }
  
  feature <- unlist(obj[[feature.unwanted]])
  
  unique.feature <- unique(feature)
  
  if (is.list(feature.time.dict)){
    feature.time <- t(as.data.frame(feature.time.dict[feature]))
    unique.feature.time <- t(as.data.frame(feature.time.dict[unique.feature]))
  } else {
    feature.time <- feature.time.dict[feature]
    unique.feature.time <- feature.time.dict[unique.feature]
  }
  
  
  # W <- 1 / (sqrt(euclidean_pdist2(x = unique.feature.time, y = feature.time)) + reg) ** pow
  W <- exp(-euclidean_pdist2(x = unique.feature.time, y = feature.time) / 2 * reg ** 2) ** (2 * pow)
  old.embedding <- obj@reductions[[reduction.use]]@cell.embeddings[, dims.use]
  #old.embedding = t(t(old.embedding) / matrixStats::colSds(old.embedding))
  
  B <- matrix(0, ncol=length(dims.use), nrow=length(dims.use))
  for (j in 1:length(unique.feature))
  {
    mask <- (feature != unique.feature[j])
    temp <- t(old.embedding[mask, ]) - colMeans(old.embedding[!mask, ])
    temp <- t(t(temp) * sqrt(W[j, mask]))
    B <- B + temp %*% t(temp) # / sum(mask)
  }
  B <- B / sum(W)
  
  cell.embeddings <- t(solve(chol(B), t(old.embedding)))
  
  if (strength > 1) {
    for (i in 2:strength) {
      cell.embeddings <- t(solve(chol(B), t(old.embedding)))
    }
  }
  
  cell.embeddings = t(t(cell.embeddings) / matrixStats::colSds(cell.embeddings))
  colnames(cell.embeddings) <- paste0(reduction.key, 1:length(dims.use))
  
  alt.reduction <- CreateDimReducObject(
    embeddings = cell.embeddings,
    key = reduction.key
  )
  obj[[reduction.name]] <- alt.reduction
  obj
}

