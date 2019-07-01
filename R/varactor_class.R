library(R6)
library(matrixStats)
library(Rtsne)
library(ggplot2)

# The core of varactor is distance.
# So there should be three fields:
# data, distance, and embedding 
# (and maybe clustering, which may just serve as a label)

# Avoid complicating the situation or playing too much with R6.
# Only use active fields when the effect is understandable (i.e. without weird side-effect).

# It should be able to hold multiple clustering/embedding results to suit its name.

# Postpone the self-assigned distance function. It's not practical to let the user to assign one after all...
# If one can, let them edit the code...

# There can be two different types of clustering methods, the ones based on distance
# and the ones not based on distance. How to solve this?

# How to define a pdist2 function? It should be a function only accept one matrix as input...

Varactor <- R6Class(
  classname = "Varactor", 
  
  public = list(
    # public fields
    
    # public methods
    initialize = function(data, labels, what = "raw", 
                          reduce = 50, verbose = TRUE, manual = FALSE){
      private$.verbose_write("Now you are creating an R6 object. 
                             If you want to make a copy for an R6 object, use 'a_copy <- this_one$clone()'.")
      
      self$verbose <- verbose

      private$.labels <- labels
      
      if (what %in% c("raw", "normalized", "reduced")) private[[paste0('.', what)]] <- data
      else stop("Bad argument what: must be 'raw', 'normalized' or 'reduced'.")
      
      if (!manual)
      {
        if (what == "raw"){
          self$normalize()
        }
        if (what == "raw" || what == "normalized"){
          self$reduce()
        }
      }
    },
    
    normalize <- function(){
      private$.verbose_write("Normalizing data...")
      private$.data$normalized <- t(t(self$.raw) / colMeans(self$.raw))
      invisible(self)
    }
    
    reduce <- function(...){
      private$.verbose_write("Calculating PCA...")
      private$.verbose_write(paste0("Dimensionality is set to ", 
                                    reduce, 
                                    ". If this is not preferred provide a value for parameter reduce."))
      private$.verbose_write("Depending on the data size, this may take a few minutes...")
      
      private$.reduced <- prcomp(t(private$.normalized), rank.=reduce, ...)$x
      invisible(self)
    }
    
    define_metric = function(name, type = "euclidean", manual = FALSE,
                               strata = NA, mahalanobis_cov = NA){
      if (type == "euclidean"){
        private$.pdist2[[name]] <- function(x){
          xx = rowSums(x * x)
          pdist <- t(matrix(rep(xx, length(xx)), ncol=length(xx)))
          pdist <- pdist + xx
          pdist <- pdist - 2 * x %*% t(x)
          pdist 
        }
      }
      
      else if (type == "mahalanobis"){
        invS <- solve(mahalanobis_cov)
        private$.pdist2[[name]] <- function(x){
          xx = rowSums(x %*% invS * x)
          pdist <- t(matrix(rep(xx, length(xx)), ncol=length(xx)))
          pdist <- pdist + xx
          pdist <- pdist - 2 * x %*% invS %*% t(x)
          pdist
        }
      }
      
      else if (type == "davidson"){
        B = matrix(0, ncol=dim(x)[2], nrow=dim(x)[2])
        if (is.na(strata)) stop("Bad argument")
        unwanted_label <- labels[[strata]]
        for (j in unique(unwanted_label))
        {
          mask <- bad_labels != j
          temp <- t(x[mask, ]) - colSums(x[!mask, ])
          B <- B + temp %*% t(temp)
        }
        B <- B / dim(x)[1]
        
        # Call this function recursively to get a mahalanobis distance metric
        # with B as its covariance matrix
        # Manual is set to TRUE in calling to avoid calculating distance twice.
        self$define_metric(name = name, type = "mahalanobis", manual = TRUE,
                           mahalanobis_cov = B)
      }
      else stop("Bad argument. Distance type not supported.")
      
      if (!manual){
        self$measure(name = name)
      }
    },
    
    measure = function(name)
    {
      temp <- private$.pdist2(private$.reduced)
      private$.distance_matrices[[name]] <- sqrt(temp - min(temp))
    },
    
    embed = function(name, type = "tsne", ...){
      if (!(name %in% names(private$.distance_matrices))){
        stop(paste0("Bad argument: distance name '", name, "' does not exist."))
      }
      
      if (type == "tsne"){
        private$.embedding[[name]] = Rtsne(private$.distance_matrices[[name]], 
                                           is_distance = T, ...)
      }
      else if (type == "umap"){
        stop("Not implemented.")
      }
      else if (type == "mds"){
        stop("Not implemented.")
      }
      else{
        stop("Bad argument: unknown embedding type")
      }
      invisible(self)
    },
    
    plot_reduced = function(label_name, dims = c(1, 2), ...){
      if (length(dims) <= 2){
        print(private$.data$reduced)
        plot(private$.data$reduced[, dims], col=factor(private$.labels[[label_name]]), ...)
        title("title")
      }
      else if (length(dims) > 2){
        pairs(private$.data$reduced[, dims], col=factor(private$.labels[[label_name]]), ...)
      }
      else stop("Bad argument")
    },
    
    
    
    plot_embedding = function(embedding_name){
      stop("Not implemented.")
    }
  ),

  private = list(
    
    # private fields
    .raw = NA,
    .normalized = NA,
    .reduced = NA,
    
    .labels = NA,
    
    .distance_matrices = list(),
    .embeddings = list(),
    
    .pdist2_func = NA,

    # verbosity
    .verbose = NA,
    .verbose_print = NA,
    .verbose_write = NA,
    
  ),
  
  active = list(
    verbose = function(value) {
      if (missing(value)) {
        return(private$.verbose)
      } 
      else {
        if (identical(value, F)){
          private$.verbose_print <- function(x, ...){invisible(x)} # supress output, but still return the same value as the normal print
          private$.verbose_write <- function(...){} # Doing nothing and return NULL (also the same as a normal writeLines)
        }
        else if (identical(value, T)){
          private$.verbose_print <- print
          private$.verbose_write <- writeLines
        }
        else stop("Bad argument: must be logical TRUE or FALSE")
        
        private$.verbose <- value
      }
    },
    
    raw = function(value){
      if (missing(value)){
        return(private$.raw)
      }
      else{
        private$.verbose_write("Note: You are changing the raw expression. 
                               This neither automatically perform preprocessing or recalculation the results. 
                               You may consider to run normalize() and reduce() manually.")
        private$.raw <- value
      }
    },
    
    normalized = function(value){
      if (missing(value)) {
        return(private$.normalized)
      } 
      else{
        private$.verbose_write("Note: You are changing the normalized expression. 
                               This neither automatically perform preprocessing or recalculation the results. 
                               You may consider to run reduce() manually. 
                               Also, discrepancy is possible between the raw data and this newly assigned normalized data.")
        private$.normalized <- value
      }
    },
    
    reduced = function(value){
      if (missing(value)) {
        return(private$.reduced)
      } 
      else{
        private$.verbose_write("Note: You are changing the normalized expression. 
                               This does not automatically recalculation the results. 
                               Also, discrepancy is possible between the raw data, normalized data and this newly assigned reduced data.")
        private$.reduced <- value
      }
    }
  )
)

