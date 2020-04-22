library(R6)
#library(matrixStats)
library(Rtsne)
library(umap)
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

# The normalization should be performed before getting the common genes...

rowSds <- function(x) sqrt(rowMeans((x - rowMeans(x)) ** 2) / dim(x)[2])

combine_labels <- function(labels, keep_sample_name = TRUE, data = NULL, data_width = NULL){
  
  if (!xor(is.null(data), is.null(data_width))){
    stop("Bad arguments: must specify one and only one of data and data_width.")
  }
  
  if (class(labels) != "list") stop("Bad argument: labels should be a list.")
  for (i in 1:length(labels)){
    if (class(labels[[i]]) != "list") stop("Bad argument: labels should contain lists.")
    for (j in 1:length(labels[[i]])){
      if (class(labels[[i]][[j]]) != "character") stop("Bad arguments: a label should be of class character")
    }
  }
  
  data_names <- names(labels)
  label_names <- unique(unlist(lapply(labels, names)))
  if (keep_sample_name){
    label_names <- c(label_names, "sample")
  }
  
  if (is.null(data_width)){
    data_width <- unlist(lapply(dataset_names, function(x) dim(data[[x]])[2]))
  }
  
  # c() if presence in all, unspecified if not provided
  names(label_names) <- label_names
  
  lapply(label_names, 
         function(x) unlist(lapply(data_names, 
                                   function(y) if (length(labels[[y]][[x]]) > 0) 
                                                 labels[[y]][[x]]
                                               else
                                                 rep(y, data_width[y])
                                  )
                           )
        )
}

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

Varactor <- R6Class(
  classname = "Varactor", 
  
  public = list(
    # public fields
    
    # public methods with side effects (all return invisible(self))
    initialize = function(data, 
                          labels, 
                          what = "raw", 
                          reduce_dim = 50, importance = "equal", # for PCA
                          verbose = TRUE, auto = FALSE # control
                          ){
      self$verbose <- verbose
      
      private$.verbose_write("Now you are creating an R6 object.")
      private$.verbose_write("If you want to make a copy for an R6 object, use 'a_copy <- this_one$clone()'.")
      
      if (what %in% c("raw", "normalized", "combined", "reduced")) private[[paste0('.', what)]] <- data
      else stop("Bad argument what: must be 'raw', 'normalized' or 'reduced'.")
      
      if (what %in% c("raw", "normalized")){
        private$.labels <- combine_labels(labels, data_width = unlist(lapply(data, function(x) dim(x)[2])))
      }
      else{
        private$.labels <- labels
      }
      
      if (auto)
      {
        if (what == "raw"){
          self$normalize()
        }
        if (what == "raw" || what == "normalized"){
          self$combine()
        }
        if (what == "raw" || what == "normalized" || what == "combined"){
          self$reduce(reduce_dim, importance)
        }
      }
    },
    
    normalize = function(cell_normalize = TRUE, 
                         log_transform = FALSE){
      private$.verbose_write("Normalizing data for each dataset...")
      private$.verbose_write("Note that the nomalization is based all the genes available in each dataset.")
      if (cell_normalize) private$.normalized <- lapply(private$.raw, function(x) t(t(x) / colMeans(x)))
      else private$.normalized <- private$.raw
      
      if (log_transform) private$.normalized <- lapply(private$.normalized, log1p)
      
      invisible(self)
    },
    
    combine = function(sd_threshold = 0.0, cv_threshold = 0.5){
      private$.verbose_write("Combine data to one single dataset...")
      
      temp <- lapply(private$.normalized, function(x) x[rowSds(x) > sd_threshold, ])
      
      temp <- lapply(temp, function(x) x[rowSds(x) / rowMeans(x) > cv_threshold, ])
      
      common_genes <- Reduce(intersect, lapply(temp, rownames))
      
      temp <- lapply(temp, function(x) x[common_genes, ])
      
      temp <- lapply(temp, function(x) (x - rowMeans(x)) / rowSds(x))
      
      for (i in temp)
        private$.verbose_write(paste(dim(i)))
      
      private$.combined <- do.call(cbind, temp)
    },
    
    reduce = function(reduce_dim = 50, importance = "equal", use_irlba = T){
      private$.verbose_write("Calculating PCA...")
      private$.verbose_write(paste0("Dimensionality is set to ", 
                                    reduce_dim, 
                                    ". If this is not preferred provide a value for parameter reduce_dim"))
      private$.verbose_write("Depending on the data size, this may take a few minutes...")
      
      gene_Sds <- rowSds(private$.combined)
      private$.verbose_write(paste0(sum(gene_Sds == 0), " constant genes are not considered in PCA."))
      
      if (importance == "equal"){
        if (!use_irlba)
          private$.reduced <- prcomp(t(private$.combined[gene_Sds > 0, ]), 
                                   rank.=reduce_dim, scale. = TRUE)$x
        else
          require(irlba)
          private$.reduced <- prcomp_irlba(x = t(private$.combined[gene_Sds > 0, ]), 
                                    n = reduce_dim, scale. = T)$x
      }
      else{
        stop("Not implemented.")
      }
      invisible(self)
    },
    
    clean = function(what = "auto"){
      if (what == "auto"){
        if (!identical(private$.normalized, NA)) private$.raw <- NA 
        if (!identical(private$.combined, NA)) private$.normalized <- NA 
        if (!identical(private$.reduced, NA)) private$.combined <- NA 
      }
      else if (what == "raw"){
        private$.raw <- NA 
      }
      else if (what == "normalized"){
        private$.normalized <- NA
      }
      else if (what == "combined"){
        private$.combined <- NA
      }
      gc()
      invisible(self)
    },
    
    define_metric = function(name, type = "euclidean", manual = TRUE,
                             strata, mahalanobis_cov, strength = 1, 
                             time_label, time_of_strata, reg=1., pow=1.){
      # Be very careful that all variables in the functions are in the local
      # environment (closure) to avoid untrackable bugs.
      if (type == "euclidean"){
        private$.metric[[name]] <- function(x){
          xx = rowSums(x * x)
          pdist <- t(matrix(rep(xx, length(xx)), ncol=length(xx)))
          pdist <- pdist + xx
          pdist <- pdist - 2 * x %*% t(x)
          pdist 
        }
      }
      
      else if (type == "canberra"){
        private$.metric[[name]] <- function(x){
          pdist <- matrix(0, ncol=dim(x)[1], nrow = dim(x)[1])
          for (i in 1:dim(x)[1]){
            for (j in i:dim(x)[1]){
              pdist[i, j] <- sum(abs(x[i, ] - x[j, ]) / (abs(x[i, ]) + abs(x[j, ])))
              pdist[j, i] <- pdist[i, j]
            }
          }
          return(pdist ** 2)
        }
      }
      
      else if (type == "mahalanobis"){
        invS <- solve(mahalanobis_cov)
        private$.metric[[name]] <- function(x){
          xx = rowSums(x %*% invS * x)
          pdist <- t(matrix(rep(xx, length(xx)), ncol=length(xx)))
          pdist <- pdist + xx
          pdist <- pdist - 2 * x %*% invS %*% t(x)
          pdist
        }
      }
      
      else if (type == "davidson_mod"){
        B <- matrix(0, ncol=dim(private$.reduced)[2], nrow=dim(private$.reduced)[2])
        if (is.na(strata)) stop("Bad argument")
        unwanted_label <- private$.labels[[strata]]
        cnt = 0;
        for (j in 1:dim(private$.reduced)[1])
        {
          for (k in 1:dim(private$.reduced)[1])
          if (unwanted_label[i] != unwanted_label[j])
          {
            cnt <- cnt + 1
            temp <- private$.reduced[i, ] - private$.reduced[j, ]
            B <- B + temp %*% t(temp)
          }
        }
        B <- B / cnt
        
        # Call this function recursively to get a mahalanobis distance metric
        # with B as its covariance matrix
        # Manual is set to TRUE in calling to avoid calculating distance twice.
        if (strength == 1)
          self$define_metric(name = name, type = "mahalanobis", manual = TRUE,
                             mahalanobis_cov = B)
        if (strength == 2)
          self$define_metric(name = name, type = "mahalanobis", manual = TRUE,
                             mahalanobis_cov = B %*% B)
      }
      
      else if (type == "davidson"){
        B <- matrix(0, ncol=dim(private$.reduced)[2], nrow=dim(private$.reduced)[2])
        if (is.na(strata)) stop("Bad argument")
        unwanted_label <- private$.labels[[strata]]
        for (j in unique(unwanted_label))
        {
          mask <- unwanted_label != j
          temp <- t(private$.reduced[mask, ]) - colMeans(private$.reduced[!mask, ])
          B <- B + temp %*% t(temp)
        }
        B <- B / dim(private$.reduced)[1]
        
        # Call this function recursively to get a mahalanobis distance metric
        # with B as its covariance matrix
        # Manual is set to TRUE in calling to avoid calculating distance twice.
        if (strength == 1)
        self$define_metric(name = name, type = "mahalanobis", manual = TRUE,
                           mahalanobis_cov = B)
        if (strength == 2)
          self$define_metric(name = name, type = "mahalanobis", manual = TRUE,
                             mahalanobis_cov = B %*% B)
      }
      
      else if (type == "continuous_davidson"){
        # naive implementation, optimization needed
        B <- matrix(0, ncol=dim(private$.reduced)[2], nrow=dim(private$.reduced)[2])
        n <- dim(private$.reduced)[1]
        
        if (is.na(strata)) stop("Bad argument")
        unwanted_label <- private$.labels[[strata]]
        
        W <- euclidean_pdist2(x = unwanted_label)
        
        for (i in 1:n){
          for (j in 1:n){
            temp <- t(private$.reduced[i, ] - private$.reduced[j, ])
            B <- B + temp %*% t(temp) * W[i, j]
          }
        }
        B <- B / sum(W)
        
        self$define_metric(name = name, type = "mahalanobis", manual = TRUE,
                           mahalanobis_cov = B)
      }
      
      else if (type == "temporal_davidson"){
        # naive implementation, optimization needed
        B <- matrix(0, ncol=dim(private$.reduced)[2], nrow=dim(private$.reduced)[2])
        n <- dim(private$.reduced)[1]
        
        if (missing(time_label)) stop("Bad argument: must provide time_label")
        unwanted_label <- time_label
        unique_unwanted_label <- unique(unwanted_label)
        W <- 1 / sqrt(euclidean_pdist2(x = unique_unwanted_label, y = unwanted_label) + 1.)
        
        for (j in 1:length(unique_unwanted_label))
        {
          mask <- unwanted_label != unique_unwanted_label[j]
          temp <- t(private$.reduced[mask, ]) - colMeans(private$.reduced[!mask, ])
          temp <- t(t(temp) * W[j, mask])
          B <- B + temp %*% t(temp)
        }
        B <- B / sum(W ** 2)
        
        self$define_metric(name = name, type = "mahalanobis", manual = TRUE,
                           mahalanobis_cov = B)
      }
      
      else if (type == "mixed_davidson"){
        # naive implementation, optimization needed
        B <- matrix(0, ncol=dim(private$.reduced)[2], nrow=dim(private$.reduced)[2])
        n <- dim(private$.reduced)[1]
        
        if (missing(time_label)) stop("Bad argument: must provide time_label for strata")
        if (missing(strata)) stop("Bad argument: must provide strata")
        unwanted_label <- private$.labels[[strata]]
        unique_unwanted_label <- unique(unwanted_label)
        
        W <- 1 / sqrt(euclidean_pdist2(x = time_of_strata, y = time_label) + reg) ** pow
        
        for (j in 1:length(unique_unwanted_label))
        {
          mask <- unwanted_label != unique_unwanted_label[j]
          temp <- t(private$.reduced[mask, ]) - colMeans(private$.reduced[!mask, ])
          temp <- t(t(temp) * W[j, mask])
          B <- B + temp %*% t(temp)
        }
        B <- B / sum(W ** 2)
        
        self$define_metric(name = name, type = "mahalanobis", manual = TRUE,
                           mahalanobis_cov = B)
      }
      
      else stop("Bad argument. Distance type not supported.")
      
      if (!manual){
        self$measure(name = name)
      }
      
      invisible(self)
    },
    
    measure = function(name){
      temp <- private$.metric[[name]](private$.reduced)
      private$.distance_matrices[[name]] <- sqrt(temp - min(temp))
      
      invisible(self)
    },
    
    embed = function(name, type = "tsne", scaling = 1., custome_func, ...){
      if (!(name %in% names(private$.distance_matrices))){
        stop(paste0("Bad argument: distance name '", name, "' does not exist."))
      }
      if (type == "customize"){
        private$.embeddings[[name]] <- func(private$.distance_matrices[[name]])
      }
      if (type == "tsne"){
        private$.embeddings[[name]] = Rtsne(private$.distance_matrices[[name]], 
                                           is_distance = T, ...)$Y
      }
      else if (type == "umap"){
        config <- umap.defaults
        config$input <- "dist"
        private$.embeddings[[name]] <- umap(private$.distance_matrices[[name]],
                                           config = config)$layout
      }
      else if (type == "umap_3d"){
        config <- umap.defaults
        config$input <- "dist"
        config$n_components <- 3
        private$.embeddings[[name]] <- umap(private$.distance_matrices[[name]],
                                            config = config)$layout
      }
      else if (type == "rbf_umap"){
        config <- umap.defaults
        config$input <- "dist"
        private$.embeddings[[name]] <- umap(sqrt(1-exp(-private$.distance_matrices[[name]] ** 2 / (2 * scaling ** 2))),
                                            config = config)$layout
      }
      else if (type == "rbf_umap_3d"){
        config <- umap.defaults
        config$input <- "dist"
        config$n_components <- 3
        private$.embeddings[[name]] <- umap(sqrt(1-exp(-private$.distance_matrices[[name]] ** 2 / (2 * scaling ** 2))),
                                            config = config)$layout
      }
      else if (type == "mds"){
        stop("Not implemented.")
      }
      else{
        stop("Bad argument: unknown embedding type")
      }
      invisible(self)
    },
    
    # public methods without side effects
    plot_reduced = function(label_name, dims = c(1, 2), ...){
      if (length(dims) <= 2){
        p <- plot(private$.reduced[, dims], col=factor(private$.labels[[label_name]]), ...)
      }
      else if (length(dims) > 2){
        p <- pairs(private$.reduced[, dims], col=factor(private$.labels[[label_name]]), ...)
      }
      else stop("Bad argument")
      
      return(p)
    },
    
    plot_embedding = function(name, label_name, ...){
      p <- plot(private$.embeddings[[name]], col=factor(private$.labels[[label_name]]), ...)
      return(p)
    },
    
    plot = function(name, what, label_name, size=1., manual = FALSE){
      if (what == 'embedding'){
        if (!manual){
          df <- data.frame(x = private$.embeddings[[name]][, 1], 
                         y = private$.embeddings[[name]][, 2],
                         l = private$.labels[[label_name]])
        
          return(ggplot(df) + geom_point(aes(x=x, y=y, color=l), size=size) + labs(color=label_name) +
                   ggtitle(paste(name, what)))
        }
        else{
          df <- data.frame(x = private$.embeddings[[name]][, 1], 
                           y = private$.embeddings[[name]][, 2])
          df[names(private$.labels)] <- private$.labels
          return(df)
        }
      }
      else{
        stop("Not implemented.")
      }
        
    }
    
  ),

  private = list(
    
    # private fields
    .raw = NA,        # a list of datasets (as matrices)
    .normalized = NA, # a list of datasets (as matrices)
    .combined = NA,   # a single combined matrix
    .reduced = NA,    # a  
    
    .raw_labels = NA,
    .labels = NA,
    
    .distance_matrices = list(),
    .embeddings = list(),
    
    .metric = list(),

    # verbosity
    .verbose = NA,
    .verbose_print = NA,
    .verbose_write = NA
    
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
        private$.verbose_write("Note: You are changing the normalized expression.")
        private$.verbose_write("This may cause discrepancy between it and the up/down stream data.")
        private$.normalized <- value
      }
    },
    
    combined = function(value){
      if (missing(value)) {
        return(private$.combined)
      } 
      else{
        private$.verbose_write("Note: You are changing the combined expression.")
        private$.verbose_write("This may cause discrepancy between it and the up/down stream data.")
        private$.combined <- value
      }
    },
    
    reduced = function(value){
      if (missing(value)) {
        return(private$.reduced)
      } 
      else{
        private$.verbose_write("Note: You are changing the reduced expression. 
                               This does not automatically recalculation the results. 
                               Also, discrepancy is possible between the raw data, normalized data and this newly assigned reduced data.")
        private$.reduced <- value
      }
    },
    
    labels = function(value){
      if (missing(value)) {
        return(private$.labels)
      } 
      else{
        private$.verbose_write("Note: You are changing the labels.")
        private$.labels <- value
      }
    },
    
    embeddings = function(value){
      if (missing(value)){
        return(private$.embeddings)
      }
      else{
        stop("Cannot edit embedding.")
      }
    },
    
    steal = function(value){
      if (missing(value)){
        return(private)
      }
      else{
        private <- value
      }
    }
  )
)

