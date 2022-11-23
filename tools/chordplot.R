# Circle Chord Plot for Tree Representation of Brain Connectomes
# This code plots circle chord plot for the Tree Representations of 
# Brain Structural Connectivity via Persistent Homology by D. Li, P. Nguyen, 
# Z. Zhang and D. Dunson
#
# Author: Phuc Nguyen, Aug. 2020

library(circlize)
library(rmatio)
library(tidyverse)
library(png)
source("tools/utils.R")
source("tools/make_tree.R")
source("tools/chordplot_adjmat.R")
source("tools/chordplot_rhlh.R")
source("tools/chordplot_lobe.R")
source("tools/chordplot_slobe.R")
source("tools/chordplot_sslobe.R")
source("tools/chordplot_tree.R")

EDGES <- read.csv("data/Desikan_Cortical_Tree_Edge_List.csv")
INDEX_MAP <- readRDS("data/Desikan_Cortical_Index_Map.Rds")
NODES <- get_tree_names(make_null_tree())
# N_REGIONS <- 68
# AMAT <- matrix(1, N_REGIONS, N_REGIONS)
MEAN_BRAIN <- read.mat("data/HCP_cortical_TensorData_desikan.mat")[["loaded_tensor"]][,,1,] %>%
  apply(c(1,2), mean)
MAT_ROI <- read.table("data/Desikan_Cortical_ROI_Order.txt", 
                  fill = TRUE, 
                  col.names = c("id", "name", "v1", "v2", "v3", "v4")) %>%
  filter(id > 1000) %>%
  separate(name, into = c("layer", "side", "desikan"), sep = "-") %>%
  mutate(roi = paste(side, desikan, sep = "-")) %>%
  dplyr::select(roi) %>%
  as.matrix() %>% as.vector()
colnames(MEAN_BRAIN) <- MAT_ROI
rownames(MEAN_BRAIN) <- MAT_ROI

# Functions
corcolfunc <- colorRampPalette(c("blue","white", "red"))

# Calculate an adjaceny matrix for the given level of the tree in id
# id: id of nodes in the tree to plot
# nodes: names of a tree vector
get_chord_mat <- function(nodes, amat, id, fun="sum") {
  parents <- c()
  factors <- c()
  # Store matrices
  Alist <- list()
  listi <- 1
  for (n in id) {
    node <- nodes[n]
    children <- get_children(node)
    l <- length(children)
    parents <- c(parents, rep(node, l))
    factors <- c(factors, children)
    if (l > 0) {
      A <- matrix(0, l, l)
      for (i in 1:(l-1)) {
        for (j in (i+1):l) {
          chi <- INDEX_MAP[[children[i]]]
          chj <- INDEX_MAP[[children[j]]]
          if (fun == "sum") {
            A[i,j] <- sum(amat[as.matrix(expand.grid(chi, chj))]) + sum(amat[as.matrix(expand.grid(chj, chi))])
          }
          if (fun == "mean") {
            A[i,j] <- mean(amat[as.matrix(expand.grid(chi, chj))]) + mean(amat[as.matrix(expand.grid(chj, chi))])
          }
          A[j,i] <- A[i,j]
        }
      }
      Alist[[listi]] <- A
      listi <- listi + 1
    }
  }
  AA <- as.matrix(Matrix::bdiag(Alist))
  colnames(AA) <- factors
  rownames(AA) <- factors
  return(AA)
}


# Calculate angle for text
xy_to_angle <- function(x, y) {
  a <- 361^(y < 0) - 1
  return(acos(x) * (180/pi) * (-1)^(y < 0) + a)
}


# Get top connections for each region from adjacency matrix
get_top_links <- function(A2) {
  A2[lower.tri(A2)] <- 0
  A2df <- data.frame(from = rep(rownames(A2), times = ncol(A2)),
                     to = rep(colnames(A2), each = nrow(A2)),
                     value = as.vector(A2),
                     stringsAsFactors = FALSE) %>%
    filter(value > 0) %>%
    arrange(desc(value))
  #A2df <- A2df[!duplicated(A2df$from) | !duplicated(A2df$to),]
  A2df <- A2df[!duplicated(A2df$from),]
  return(A2df)
}

# Manually create representative links in the lobe level here for visual purpose.
get_lobe_links <- function() {
  data.frame(
    from = c(L3_BAC_FRO_RIG[1],
             L3_BAC_FRO_RIG[2],
             L3_BAC_FRO_RIG[4],
             L3_BAC_FRO_RIG[5],
             L3_FRO_BAC_LEF[6],
             L3_FRO_BAC_LEF[5],
             L3_FRO_BAC_LEF[3],
             L3_FRO_BAC_LEF[2]),
    to = c(L3_BAC_FRO_RIG[3],
           L3_BAC_FRO_RIG[4],
           L3_BAC_FRO_RIG[5],
           L3_BAC_FRO_RIG[6],
           L3_FRO_BAC_LEF[4],
           L3_FRO_BAC_LEF[3],
           L3_FRO_BAC_LEF[2],
           L3_FRO_BAC_LEF[1]),
    value = 10
  )
}

# Get top n links with largest values
get_sig_links <- function(A2, topn=NULL) {
  A2[lower.tri(A2)] <- 0
  A2df <- data.frame(from = rep(rownames(A2), times = ncol(A2)),
                     to = rep(colnames(A2), each = nrow(A2)),
                     value = as.vector(A2),
                     stringsAsFactors = FALSE) %>%
    filter(value != 0)
  if (!is.null(topn)) {
    A2df <- A2df %>%
      arrange(desc(abs(value))) %>%
      head(topn)
  }
  return(A2df)
}


# Get ribbon size
# Size for L3_BAC_FRO
LOBE_RIBSIZE_MAP <- c("rh-insula"   = 6    ,
                      "rh-temporallobe" = 24 ,
                      "rh-cingulate"  = 12  ,
                      "rh-occipitallobe" = 14,
                      "rh-parietallobe" = 14,
                      "rh-frontallobe"  = 30)
LOBE_RIBSIZE_BAC_FRO <- LOBE_RIBSIZE_MAP[L3_BAC_FRO_RIG]

LOBE_RIBSIZE <- function() {
  #return(c(c(10, 6, 4, 10, 28, 42), rev( c(10, 6, 4, 10, 28, 42))))
  return(c(LOBE_RIBSIZE_BAC_FRO, rev(LOBE_RIBSIZE_BAC_FRO)))
}

SUBLOBE_RIBSIZE <- function(gap.after=1.5) {
  midid <- (length(L4)/2)
  ribn <- sapply(L3, function(n) max(1, length(get_children(n))))
  ribsize <- (LOBE_RIBSIZE() - gap.after*0.5*(ribn - 1)) / (ribn)
  ribsize <- unlist(sapply(1:length(ribn), 
                           function(i) rep(ribsize[i], ribn[i])))
  return(ribsize)
}

SSLOBE_RIBSIZE <- function(gap.after=1.5) {
  midid <- (length(L5)/2)
  ribn <- sapply(L4, function(n) max(1, length(get_children(n))))
  gap.after <- 1
  ribsize <- (SUBLOBE_RIBSIZE() - gap.after*0.75*(ribn - 1)) /(ribn)
  ribsize <- unlist(sapply(1:length(ribn), 
                           function(i) rep(ribsize[i], max(ribn[i], 1))))
  return(ribsize)
}

LEAVES_RIBSIZE <- function() {
  midid <- (length(L6)/2)
  ribn <- sapply(L5, function(n) length(get_children(n)))
  ribsize <- SSLOBE_RIBSIZE() /
    (ribn + 1)
  ribsize <- unlist(sapply(1:length(ribn), 
                           function(i) rep(ribsize[i], max(ribn[i], 1))))
  return(ribsize)
}

L2 <- c("rh", "lh")
L3_BAC_FRO_RIG <- c("rh-insula"       ,
                    "rh-temporallobe" ,
                    "rh-cingulate"    ,
                    "rh-occipitallobe",
                    "rh-parietallobe" ,
                    "rh-frontallobe"  )
L3_FRO_BAC_LEF <- c("lh-frontallobe"  ,
                    "lh-parietallobe" ,
                    "lh-occipitallobe",
                    "lh-cingulate"    ,
                    "lh-temporallobe" ,
                    "lh-insula")
L3 <- c(L3_BAC_FRO_RIG, L3_FRO_BAC_LEF)
L4 <- unlist(sapply(L3, function(n) {
  if (length(grep("rh", n)) > 0) {ch <- rev(get_children(n))}
  else {ch <- get_children(n)}
  if (length(ch)==0) {ch <- n}
  return(ch)
  }))
L5 <- unlist(sapply(L4, function(n) {
  if (length(grep("rh", n)) > 0) {ch <- rev(get_children(n))}
  else {ch <- get_children(n)}
  if (length(ch)==0) {ch <- n}
  return(ch)
}))
L6 <- unlist(sapply(L5, function(n) {
  if (length(grep("rh", n)) > 0) {ch <- rev(get_children(n))}
  else {ch <- get_children(n)}
  if (length(ch)==0) {ch <- n}
  return(ch)
}))
L2ID <- c(2, 13) # id for rh, lh

REGION_IMG_MAP <- list("rh-insula" = "dkk-right-insula.png",
                       "rh-temporallobe" = "dkk-right-temporal.png" ,
                       "rh-cingulate" = "dkk-right-cingulate.png"   ,
                       "rh-occipitallobe" = "dkk-right-occipital.png",
                       "rh-parietallobe" = "dkk-right-parietal.png",
                       "rh-frontallobe" = "dkk-right-frontal.png",
                       "lh-frontallobe" = "dkk-left-frontal.png" ,
                       "lh-parietallobe" = "dkk-left-parietal.png",
                       "lh-occipitallobe" = "dkk-left-occipital.png",
                       "lh-cingulate" = "dkk-left-cingulate.png"   ,
                       "lh-temporallobe" = "dkk-left-temporal.png",
                       "lh-insula" = "dkk-left-insula.png"
                       )