
# Methods ---------------------------------------------------------------------

## k means ====================================================================
#' Get k-means Clustering
#'
#' This is just a wrapper funciton for the kmeans function
#'
#' @title Get k-means Clustering
#' @aliases
#'
#' @author Daniel Kick (\email{daniel.r.kick@@protonmail.com})
#'
#' @keywords Clustering
#'
#' @export
#'
#' @return clustering assignments
#'
#' @examples
#'
get_kmeans_clustering <- function(input.df = mrna_raw[, -1],
                                  target.nclusters = 11) {
  kmeans.m <- kmeans(input.df, centers = target.nclusters)
  return(kmeans.m$cluster)
}

## H cluster ==================================================================
#' Get Hierarchical Clustering
#'
#' This is just a wrapper funciton for the pvclust function. Uses dendextend::cutree to get the target clusters.
#' See also https://cran.r-project.org/web/packages/dendextend/vignettes/Cluster_Analysis.html it's a great reference!
#'
#' @title Get Hierarchical Clustering
#' @aliases
#'
#' @author Daniel Kick (\email{daniel.r.kick@@protonmail.com})
#'
#' @keywords Clustering
#'
#' @export
#'
#' @return clustering assignments
#'
#' @examples
#'
get_hierarchical_clustering <- function(input.df = mrna_target[, -1],
                                        target.nclusters = 12,
                                        use.method.dist = "cor",
                                        use.method.hclust = "ward.D",
                                        use.nboot = 10) {
  temp <- input.df

  temp.clust <- pvclust(t(temp),
    method.dist = use.method.dist,
    method.hclust = use.method.hclust,
    nboot = use.nboot
  )
  temp.clust <- as.dendrogram(temp.clust$hclust)

  assignment <- cutree(temp.clust, k = target.nclusters)[order.dendrogram(temp.clust)]

  return(assignment)
}
# hclust_methods <- c("ward.D", "single", "complete", "average", "mcquitty", "median", "centroid", "ward.D2")
# dist_methods <- c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", #"correlation", "uncentered")


## SNN-Cliq ===================================================================
#' Get SNN-Cliq Clustering
#'
#' This is just a wrapper funciton for the SNN function. Note that at the time of writing this has not been set up to work on linux/macos.
#'
#' @title Get SNN Clustering
#' @aliases
#'
#' @author Daniel Kick (\email{daniel.r.kick@@protonmail.com})
#'
#' @keywords Clustering
#'
#' @export
#'
#' @return clustering assignments
#'
#' @examples
#'

# 2 tuning parameters (k, distance)
get_SNN_clustering <- function(input.df = mrna_raw[, -1],
                               use.k = 3,
                               use.distance = "euclidean") {
  temp <- input.df
  SNN(temp,
    outfile = paste0(getwd(), "/inst/extdata/output_files/temp_edges.txt"),
    k = use.k,
    distance = use.distance
  )

  if (Sys.info()["sysname"] == "Windows") {
    # in CMD run:
    # Py\Cliq.py -i inst\extdata\output_files\temp_edges.txt -o inst\extdata\output_files\temp_clust.txt
    shell("Py\\Cliq.py -i inst\\extdata\\output_files\\temp_edges.txt -o inst\\extdata\\output_files\\temp_clust.txt")
  } else if (Sys.info()["sysname"] == "Darwin") {
    warning("get_SNN_clustering() hasn't been tested on macos!")
    # TODO edit this to work on unix systems
  } else {
    warning("get_SNN_clustering() hasn't been tested on linux!")
    # TODO edit this to work on unix systems
  }

  assignment <- read.table(paste0(getwd(), "/inst/extdata/output_files/temp_clust.txt"),
    header = FALSE, sep = "", dec = "."
  )
  return(assignment)
}

# Evaluation ------------------------------------------------------------------
## Evaluate clustering performance ============================================
#' Get Cluster comparisons
#'
#' This is function uses NMF::purity, clues::adjustedRand, and BiocGenerics::as.vector to produce several common coherence measures.
#'
#' @title Evaluate Cluster Performance
#' @aliases
#'
#' @author Daniel Kick (\email{daniel.r.kick@@protonmail.com})
#'
#' @keywords Clustering, Evaluate
#'
#' @export
#'
#' @return Coherence metrics
#'
#' @examples
#'
get_cluster_comparisons <- function(reference.clustering = mrna_raw$Cell,
                                    generated.clustering) {
  if (length(reference.clustering) != length(generated.clustering)) {
    warning("Input vectors are not of the same length!")
  } else {
    # reference.clustering = mrna_raw$Cell
    # generated.clustering = kmeans.m$cluster
    output <- array(0, dim = 6)

    reference.clustering <- as.numeric(reference.clustering) %>% as.factor()
    generated.clustering <- as.numeric(generated.clustering) %>% as.factor()

    output <- NMF::purity(reference.clustering, generated.clustering)
    names(output) <- "Purity"
    # Get a lot of concurrance measures
    # https://davetang.org/muse/2017/09/21/adjusted-rand-index/
    output <- c(
      output,
      clues::adjustedRand(
        BiocGenerics::as.vector(reference.clustering),
        BiocGenerics::as.vector(generated.clustering)
      )
    )
    return(output)
  }
}

## sweep over free clustering parameters ======================================
#' Evaluate Cluster Parameters
#'
#' This uses previously defined utility functions to easily run many clustering methods and evaluate them against a known reference
#'
#' @title Evaluate Cluster Parameters
#' @aliases
#'
#' @author Daniel Kick (\email{daniel.r.kick@@protonmail.com})
#'
#' @keywords Clustering, Evaluation
#'
#' @export
#'
#' @return clustering scores for methods and parameter sets tested.
#'
#' @examples
#'
evaluate_clust_params <- function(use.input.df = mrna_raw[,-1], # mrna_cell #mrna_target #mrna_raw
                                  use.k.param.kmeans = c(11),
                                  use.k.param.hclust = c(11),
                                  use.k.param.snncliq = c(3:9), # number of neighbors to consider
                                  use.reference.ids = mrna_raw$Cell,
                                  use.dist.methods.hclust = c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "correlation", "uncentered"),
                                  use.methods.hclust = c("ward.D", "single", "complete", "average", "mcquitty", "median", "centroid", "ward.D2"),
                                  use.nboot.hclust = 1,
                                  use.dist.methods.snncliq =  c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"),
                                  reduce.clusterings = FALSE

) {
  # Generate Clusterings --------------------------------------------------------
  ## K Means Clustering =========================================================

  # use.k.param <- c(11)
  temp <- map(use.k.param.kmeans, function(iter.k.param) {
    get_kmeans_clustering(
      input.df = use.input.df,
      target.nclusters = iter.k.param
    )
  })

  param.combinations <- paste("K", as.character(
    rep(use.k.param.kmeans, each = length(use.k.param.kmeans))
  ), sep = ".")

  K.clusters <- do.call(cbind.data.frame, temp)
  names(K.clusters) <- param.combinations

  ## Hierarchical Clustering ====================================================
  #use.dist.methods <- use.dist.methods.hclust #c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "correlation", "uncentered")
  #use.hclust.methods <- use.hclust.methods.hclust #c("ward.D", "single", "complete", "average", "mcquitty", "median", "centroid", "ward.D2")

  temp <- map(use.dist.methods.hclust, function(iter.dist) {
    map(use.methods.hclust, function(iter.hclust) {
      get_hierarchical_clustering(
        input.df = use.input.df,
        target.nclusters = use.k.param.hclust,
        use.method.dist = iter.dist, #
        use.method.hclust = iter.hclust, #
        use.nboot = use.hclust.nboot
      )
    })
  })

  param.combinations <- paste("H", as.character(
    rep(use.dist.methods.hclust, each = length(use.methods.hclust))
  ), as.character(
    rep(use.methods.hclust, times = length(use.dist.methods.hclust))
  ), sep = ".")

  H.clusters <- do.call(cbind.data.frame, temp)
  names(H.clusters) <- param.combinations
  ## SNN Clustering =============================================================
  #c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski") # Note this can't use cor and uncentered like Hclust can.

  temp <- map(use.dist.methods.snncliq, function(iter.dist) {
    map(use.k.param.snncliq, function(iter.k.param) {
      print(paste(as.character(iter.dist), as.character(iter.k.param)))

      get_SNN_clustering(
        input.df = use.input.df,
        use.k = iter.k.param, # number of neigbhors should go up to 9 for scPCRc
        use.distance = iter.dist
      )
    })
  })

  param.combinations <- paste("SNN", as.character(
    rep(use.dist.methods.snncliq, each = length(use.k.param.snncliq))
  ), as.character(
    rep(use.k.param.snncliq, times = length(use.dist.methods.snncliq))
  ), sep = ".")

  SNN.clusters <- do.call(cbind.data.frame, temp)
  names(SNN.clusters) <- param.combinations

  # Evaluate and summarize clusters ---------------------------------------------
  ## Merge Data Frames ==========================================================
  #use.algos <- c("kmeans", "hclust", "snncliq")
  #if ("kmeans" %in% use.algos & "hclust" %in% use.algos)

  all.clusters <- cbind(K.clusters, H.clusters) %>% cbind(SNN.clusters)

  if (reduce.clusterings == TRUE){
    # some clusterings fail and assing all values zero
    # Some have a TON of clusters
    num.clusters <- map(all.clusters, function(x) {
      print(max(x))
    }) %>% unlist() # get number of clusters in each group
    real.num <- use.reference.ids %>% as.numeric() %>% max()
    selection.vector <- (num.clusters < (2 * (real.num)) & num.clusters > (floor(real.num / 2)) & num.clusters > 1) # don't consider anything above 2X the real number.
    all.clusters <- all.clusters[, selection.vector]
  }

  # added for comparison; positive control
  all.clusters <- cbind(Cell.Type = as.numeric(use.reference.ids), all.clusters[])

  ## Score Clusterings ==========================================================
  temp <- map(all.clusters, function(iter.param) {
    get_cluster_comparisons(reference.clustering = use.reference.ids, generated.clustering = iter.param)
  })

  all.clusters.scores <- do.call(cbind.data.frame, temp)
  all.clusters.scores <- t(all.clusters.scores) %>% as.data.frame()
  all.clusters.scores$Clustering <- row.names(all.clusters.scores)
  #all.clusters.scores.long <- gather(all.clusters.scores, Metric, Value, 1:6)

  return(all.clusters.scores)
}

