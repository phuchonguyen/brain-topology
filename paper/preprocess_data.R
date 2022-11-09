# Libraries
library(rmatio)
library(R.matlab)
library(tidyverse)
source("tools/make_tree.R")

# Indices for type of connections in connectoms
FIBER_COUNT <- 1

# Load data
brain <- read.mat("data/HCP_cortical_TensorData_desikan.mat")  # R.matlab is too slow with this
traits <- readMat("data/HCP_Covariates.mat")  # R.matlab works with this but rmatio doesn't

# Convert connectomes into trees
# brain$mode3[["count"]] is the number of connecting fibers
feat <- FIBER_COUNT
ctm <- apply(brain$loaded_tensor[,,feat,], 3, make_tree)
ctm <- t(ctm)
rownames(ctm) <- brain$all_id
saveRDS(ctm, file = "data/Desikan_Cortical_Tree_Fiber_Weights.Rds")
writeMat("data/trees.mat", trees = ctm)


# Align traits observations with connectomes
traits_df <- data.frame(id = traits$allsubject.id,
                   as.array(t(traits$cog.measure), dim = dim(t(traits$cog.measure))),
                   as.array(t(traits$confond.variable), dim = dim(t(traits$confond.variable))),
                   handness = as.matrix(traits$handness))
brain_df <- data.frame(id = brain$all_id)
traits_alg <- left_join(brain_df, traits_df, by = "id")
saveRDS(traits_alg[,2:(nrow(traits$cog.measure) + 1)],
        file = "data/CogMeasure_Aligned.Rds")
saveRDS(traits_alg[,(nrow(traits$cog.measure) + 2):(nrow(traits$cog.measure) + 4)],
        file = "data/Confond_Aligned.Rds")
saveRDS(traits_alg[,1],
        file = "data/ID_Aligned.Rds")
saveRDS(traits_alg["handness"],
        file = "data/Handness_Aligned.Rds")
