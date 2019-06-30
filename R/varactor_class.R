library(R6)
library(matrixStats)
library(Rtsne)

Varactor <- R6Class(
  classname = "Varactor", 
  
  public = list(
    initialize = function(dataset, label, what = "raw", dim = 50, verbose = T){
      self$verbose <- verbose
      self$flag <- TRUE
      self$label <- label
      if (what %in% c("raw", "normalized", "reduced")) self$data[[what]] <- dataset
      else stop("Bad argument what: must be 'raw', 'normalized' or 'reduced'")
      
      if (what == "raw"){
        self$data$log_normalized <- t(t(dataset) / colMeans(dataset))
      }
      if (what == "raw" || what == "normalized"){
        self$verbose_print("Calculating PCA. This may take a few minutes...")
        self$data$reduced <- prcomp(t(self$data$normalized), rank.=dim)$x
      }
    },
    use_pdist2_func = function(type = "euclidean", ...){
      self$pdist2 <- function(x)
      {
        xx = rowSums(x * x)
        pdist <- t(matrix(rep(xx, length(xx)), ncol=length(xx)))
        pdist <- pdist + xx
        pdist <- pdist - 2 * x %*% t(x)
        pdist
      }
    }
    set_embed_func = function(type){
      
      if (type == "tsne"){self$embed_func = function(X) Rtsne(X, is_distance = T)}
      else stop("Bad argument.")
    }
    
    plot_embedding = function(type){
      
    }
  ),
  
  private = list(
    embedding_is_valid = TRUE,
    
    data = list(raw = NA, normalized = NA, reduced = NA),
    label = NA,
    
    embedding = NA,
    
    verbose_print = NA,
    pdist2_func = NA,
    embed_func = NA
  ),
  
  active = list(
    verbose = function(value) {
      if (missing(value)) {
        return(identical(self$verbose_print, print))  
      } 
      else {
        if (identical(value, T)){
          self$verbose_print <- function(x, ...){invisible(x)}
        }
        else if (identical(value, F)){
          self$verbose_print <- print
        }
        else stop("Bad argument: must be logical TRUE or FALSE")
      }
    }
  )
)

