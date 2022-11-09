# Libraries
library(dplyr)
source("tools/utils.R")


# Load edge list for tree representation
#---------------------------------------------------------------------
# HC Connectomes based on Desikan protocol
#---------------------------------------------------------------------
EDGES <- read.csv("data/Desikan_Cortical_Tree_Edge_List.csv")
NAME <- unique(c(as.character(EDGES$from), as.character(EDGES$to)))
# Load index map, see make_index_map.R
INDEX_MAP <- readRDS("data/Desikan_Cortical_Index_Map.Rds")


# Creat a tree. 
# Each node's weight is the number of connections between its children
# A is a triangular adjacency matrix

make_tree <- function(A) {
  # Place-holder for tree
  atree <- c()
  for (node in NAME) {
    children <- get_children(node)
    l <- length(children)
    if (l == 0) { # When node is a leaf
      k <- INDEX_MAP[[node]]
      atree <- c(atree, setNames(A[k, k], node))
    } 
    else {
      s <- 0
      for (i in 1:(l - 1)) {
        for (j in (i+1):l) {
          chi <- INDEX_MAP[[children[i]]]
          chj <- INDEX_MAP[[children[j]]]
          s <- s + sum(A[chi, chj]) + sum(A[chj, chi]) # A is an triangular matrix
        }
      }
      atree <- c(atree, setNames(s, node))  
    }
  }
  
  return(atree)
}

# Helper functions for this tree
get_tree_names <- function(atree) {
  return(names(atree))
}

get_tree_weights <- function(atree) {
  return(as.vector(atree))
}

set_tree_weights <- function(atree, w) {
  names <- get_tree_names(atree)
  return(replace(atree, names, w))
}

make_null_tree <- function() {
  atree <- rep(0, length(NAME))
  names(atree) <- NAME
  return(atree)
}