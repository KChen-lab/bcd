
mahalanobis_pdist2 <- function(x, y = NULL, B = NULL, t=1)
{
  x <- t(t(x))
  if (is.null(y)) y = x
  if (is.null(B)) B = diag(dim(x)[2])
  invB <- solve(B) %^% t
  xx = rowSums(x %*% invB * x)
  yy = rowSums(y %*% invB * y)
  pdist <- t(matrix(rep(yy, length(xx)), ncol=length(yy)))
  pdist <- pdist + xx
  pdist <- pdist - 2 * x %*% invB %*% t(y)
  pdist
}

bad_corr <- function(x, bad_labels)
{
  B = matrix(0, ncol=dim(x)[2], nrow=dim(x)[2])
  for (j in unique(bad_labels))
  {
    mask <- bad_labels != j
    temp = t(x[mask, ]) - colSums(x[!mask, ])
    B <- B + temp %*% t(temp)
  }
  B <- B / dim(x)[1]
  B
}

varactor <- function(x, bad_labels, strength = 1, pca_dim = 50, pca = TRUE)
{
  if (pca)
  {
    pca_res <- prcomp(t(expr_matrix), rank.=pca_dim)
    x <- pca_res$x
  }
  B <- bad_corr(x, bad_labels)
  pdist <- sqrt(mahalanobis_pdist2(x, B=B, t=strength))
  pdist
}