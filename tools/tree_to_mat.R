# Library
library(dplyr)
source("tools/utils.R")

# Load index map, see make_index_map.R
INDEX_MAP <- readRDS("data/Desikan_Cortical_Index_Map.Rds")
EDGES <- read.csv("data/Desikan_Cortical_Tree_Edge_List.csv")
LEAVES <- unique(as.vector(unlist(INDEX_MAP)))
D <- length(LEAVES)

# Helper

# Convert a tree into an adjacency matrix
# lower: Default F returns a symmetric matrix, T returns a lower triangular matrix
make_mat <- function(atree, lower = F) {
  nm <- get_tree_names(atree)
  wt <- get_tree_weights(atree) 
  n <- length(nm)
  id_names <- c()
  A <- matrix(0, D, D)
  for (i in 1:n) {
    node <- nm[i]
    weight <- wt[i]
    children <- get_children(node)
    l <- length(children)
    if (l == 0) { # When node is a leaf
      k <- INDEX_MAP[[node]]
      A[k,k] <- weight
      id_names <- c(id_names, node)
    } 
    else {
      for (r in 1:(l - 1)) {
        for (s in (r+1):l) {
          chr <- INDEX_MAP[[children[r]]]
          chs <- INDEX_MAP[[children[s]]]
          A[chr, chs] <- weight
          A[chs, chr] <- weight
        }
      }
    }
  }
  if (lower) {
    A <- A[lower.tri(A, diag = T)]
  }
  colnames(A) <- id_names
  rownames(A) <- id_names
  return(A)
}

