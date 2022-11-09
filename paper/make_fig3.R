library(dplyr)
library(viridis)
source("tools/chordplot.R")
source("tools/make_tree.R")
source("tools/tree_to_mat.R")

# Load traits
cog <- readRDS("data/CogMeasure_Aligned.Rds")

# Load adjacency matrices
nets_mat <- readMat("data/matrices.mat")[[1]]

# Given value [0,1], return an RGB value for a colour on the given ramp
map_colors<-colorRamp(viridis(50))

# Normalize
zscores <- function(y){ qnorm( rank(y)/( length(y)+1 ) ) }
# Plot matrix
plot_matrix <- function(A, legend.name="normalized\ndifference\nfrom the mean", title="") {
  p <- A %>% 
    reshape2::melt() %>% 
    ggplot(aes(Var1, Var2, fill = zscores(value))) + 
    geom_tile() + 
    scale_fill_viridis(name = legend.name) + 
    theme(axis.title = element_blank(), 
          axis.text.x = element_text(angle = 90, size = 5), 
          axis.text.y = element_text(size = 5)) + 
    labs(title = title) + 
    theme(axis.text.x = element_blank(), 
          axis.text.y = element_blank(), 
          axis.ticks = element_blank(),
          rect = element_rect(fill = "transparent") # all rectangles
          ) + 
    coord_equal()
  return(p)
}

plot_b_tree <- function(A_tree, title="") {
  map_colors <- colorRamp(viridis(50))
  A_tree_scaled <- (A_tree - min(A_tree))/(max(A_tree) - min(A_tree))
  colorpal <- map_colors(A_tree_scaled)
  colorpal <- apply(colorpal, 1, function(x) rgb(x[1], x[2], x[3], alpha=225, maxColorValue=255))
  names(colorpal) <- names(A_tree)
  chord_tree(ribcol = colorpal, legend = T, tree = A_tree, 
             title = title,
             palette = "viridis", legend.x = 1.05, 
             img.labeled = T)
  # chord_tree(ribcol = colorpal, legend = T, tree = A_tree, 
  #            title = title,
  #            palette = "viridis", legend.x = 1.2, label.choice = c(F,T,F,F))
}

transparent_bg <- rgb(1, 1, 1, alpha=0)

# Mean brain
name <- "mean"
null_tree <- make_null_tree()
w <- get_tree_weights(null_tree)
a_tree <- set_tree_weights(null_tree, w)
a_mat <- make_mat(a_tree)
a_mat <- a_mat + t(a_mat) + (MEAN_BRAIN + t(MEAN_BRAIN))
plot_matrix(a_mat)
ggsave(paste0("figures/paper/fig3/", name, "-mat.svg"), 
       width = 7, height = 5.5, bg = "transparent")
b_tree <- make_tree(a_mat)
svglite::svglite(paste0("figures/paper/fig3/", name, "-tree.svg"), 
                 width = 7, height = 5.5,
                 bg=transparent_bg)
plot_b_tree(b_tree, title=name)
dev.off()

# Temporal lobe
name <- "temporal"
null_tree <- make_null_tree()
w <- get_tree_weights(null_tree)
#i <- c(6,7,17,18) # temporal lobe
#w[i] <- 30000
#i <- c(which(grepl("temporal", names(null_tree)))) # frontal
#w[i] <- 5000
a_tree <- set_tree_weights(null_tree, w)
a_mat <- make_mat(a_tree)
a_mat <- a_mat + t(a_mat) + (MEAN_BRAIN + t(MEAN_BRAIN))
plot_matrix(a_mat)
ggsave(paste0("figures/paper/fig3/", name, "-mat.svg"), 
       width = 7, height = 5.5, bg = "transparent")
b_tree <- make_tree(a_mat)
svglite::svglite(paste0("figures/paper/fig3/", name, "-tree.svg"), 
                 width = 7, height = 5.5,
                 bg=transparent_bg)
plot_b_tree(b_tree, title=name)
dev.off()

# Visualize some trees and corresponding adj matrices
#j <- c(which(grepl("frontal", names(null_tree)))) # frontal
#w[j] <- 60000
#w[which(grepl("frontallobe", names(null_tree)))] <- 3000
#j <- c(which(grepl("occipital", names(null_tree))))
#w[j] <- 15000
# w[which(names(null_tree)=="ctx")] <- 50
#w[which(names(null_tree) %in% c("lh", "rh"))] <- 50
a_tree <- set_tree_weights(null_tree, w)
a_mat <- make_mat(a_tree)
a_mat <- a_mat + t(a_mat) + (MEAN_BRAIN + t(MEAN_BRAIN))
plot_matrix(a_mat)
ggsave(paste0("figures/paper/fig3/frontal-mat.svg"), 
       width = 7, height = 5.5, bg = "transparent")
b_tree <- make_tree(a_mat)
svglite::svglite(paste0("figures/paper/fig3/frontal-tree.svg"), 
                 width = 7, height = 5.5,
                 bg=transparent_bg)
plot_b_tree(zscores(b_tree), title="frontal")
dev.off()


# Looking at the top and bottom quantile of some traits (images in archive)
map_colors <- colorRamp(c("#0072B2",
                          "#D3D3D3", #gray #"#FFFFFF", # white
                          "#D55E00"))
# Trait id
trait_id <- 25
trait_name <- "spatial"
trait_id <- 1; trait_name <- "fluid"
trait_id <- 7; trait_name <- "picvocab"
trait_id <- 5; trait_name <- "oralreading"

# Brain of the person who scores highest in Oral Reading Recognition
# A_tri <- MEAN_BRAIN - nets_mat[,, which(cog[,5] == max(cog[,5]))]
A_bot <- apply(nets_mat[,, which(cog[,trait_id] >= quantile(cog[,trait_id], 0.1, na.rm=T))], 
               c(1,2), mean)
A_tri <- apply(nets_mat[,, which(cog[,trait_id] >= quantile(cog[,trait_id], 0.9, na.rm=T))], 
                            c(1,2), mean) - 
    A_bot
    #MEAN_BRAIN
# A_tri <- MEAN_BRAIN - nets_mat[,, which(cog[,5] == min(cog[,5]))]
# A_tri <- apply(nets_mat[,, which(cog[,trait_id] <= quantile(cog[,trait_id], 0.1, na.rm=T))], 
#                             c(1,2), mean) - MEAN_BRAIN
A <- (A_tri + t(A_tri)) / (A_bot + t(A_bot) + 1)
colnames(A) <- MAT_ROI
rownames(A) <- MAT_ROI
# axis_ticks <- sapply(MAT_ROI, get_parent)
# axis_ticks[which(axis_ticks=="rh")] <- "rh-insula"
# axis_ticks[which(axis_ticks=="lh")] <- "lh-insula"
# axis_ticks <- sapply(axis_ticks, function(x) str_replace(x, "medialaspect|lateralaspect", "temporallobe"))
# axis_ticks <- sapply(axis_ticks, function(x) str_replace(x, "middle|inferior|orbito", ""))
# axis_ticks[duplicated(axis_ticks)] <- ""

#(MEAN_BRAIN + t(MEAN_BRAIN))
A %>% 
  reshape2::melt() %>% 
  ggplot(aes(Var1, Var2, fill = value)) + 
    geom_tile() + 
    scale_fill_gradient2(low = "#0072B2", mid = "#D3D3D3", 
                         high = "#D55E00", midpoint = 0,
                         name = "") +
    #scale_fill_viridis(name = "difference\nfrom the mean") + 
    theme(axis.title = element_blank(), 
          axis.text.x = element_text(angle = 90, size = 5), 
          axis.text.y = element_text(size = 5)) + 
  labs(title = paste("Top 10%", trait_name)) + 
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks = element_blank(),
        rect = element_rect(fill = "transparent") # all rectangles
        ) + 
  coord_equal() 
  # guides(fill = F)
ggsave(paste0("figures/paper/fig3/", trait_name, "-mat.svg"), 
       width = 7, height = 5.5, bg = "transparent")

mean_tree <- make_tree(MEAN_BRAIN)[1:23]
# top_tree <- make_tree(nets_mat[,, which(cog[,5] == max(cog[,5]))])[1:23]  # only first 23 nodes have non-zero values
# bot_tree <- make_tree(nets_mat[,, which(cog[,5] == min(cog[,5]))])[1:23]
top_tree <- make_tree(apply(nets_mat[,, which(cog[,trait_id] >= quantile(cog[,trait_id], 0.9, na.rm=T))], 
                            c(1,2), 
                            mean, na.rm=T))[1:23]
bot_tree <- make_tree(apply(nets_mat[,, which(cog[,trait_id] <= quantile(cog[,trait_id], 0.1, na.rm=T))], 
                            c(1,2), 
                            mean, na.rm=T))[1:23]

A_tree <- (top_tree - bot_tree)/bot_tree #(bot_tree - mean_tree)/mean_tree
#A_tree <- bot_tree - mean_tree #(bot_tree - mean_tree)/mean_tree
# A_tree_scaled <- (A_tree - min(A_tree))/(max(A_tree) - min(A_tree))
# colorpal <- map_colors(A_tree_scaled)
# colorpal <- apply(colorpal, 1, function(x) rgb(x[1], x[2], x[3], maxColorValue=255))
A_tree_scaled <- my_scale(A_tree, type="polar")
colorpal <- map_colors(A_tree_scaled)
colorpal <- apply(colorpal, 1, function(x) rgb(x[1], x[2], x[3],alpha=150,maxColorValue=255))
names(colorpal) <- names(A_tree)
svglite::svglite(paste0("figures/paper/fig3/",trait_name,"-tree.svg"), 
                 width = 7, height = 5.5,
                 bg=transparent_bg)
chord_tree(ribcol = colorpal, legend = T, tree = A_tree, 
           title = paste("Top 10%", trait_name),
           palette = "redwhiteblue", legend.x = 1.05, 
           img.labeled = T)
dev.off()
  # chord_tree(ribcol = colorpal, legend = T, tree = A_tree, 
#            title = paste("Bottom 10%", trait_name),
#            palette = "viridis", legend.x = 1.2, label.choice = c(F,T,F,F))
# Not very useful here. People who score high will have more connections in general.
# The nodes higher up in the tree will have more connections. 
# So visualizing raw connections is not useful. 
# The percentage change in connection doesn't show clear patterns
