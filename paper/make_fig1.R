# Make figure to explain what the tree is

library(brainconn)
library(ggplot2)
library(ggseg)
source("tools/make_tree.R")
source("tools/tree_to_mat.R")
source("tools/chordplot.R")
nodes <- get_tree_names(make_null_tree())
colournames <- unique(c("ctx", L2, L3, L4))
colournames <- colournames[which(colournames %in% nodes[1:23])]
colournames <- c(colournames, nodes[which(!nodes %in% colournames)])
colournames <- str_replace_all(colournames, "rh", "llh")
colournames <- str_replace_all(colournames, "^lh", "rh")
colournames <- str_replace_all(colournames, "llh", "lh")
colourpal <- c(viridis::viridis(n = 23, option = "D"), rep("blue", 91-23))
#colourpal <- c(viridis::plasma(n = 23), rep("gray", 91-23))
names(colourpal) <- colournames
# Make chord plot
# Superimpose these plots on each other
chord_rhlh(ribcol=colourpal[1], labeled = F)
par(new=TRUE)
chord_lobe(ribcol=colourpal, labeled = F,
           h.ratio=0.4, img.labeled = T,
           gap.after = 1.5)
par(new=TRUE)
chord_sublobe(ribcol=colourpal, labeled = F,
              h.ratio=0.65, 
              gap.after = 1.5)
par(new=TRUE)
chord_sslobe(ribcol=colourpal, labeled = F,
             h.ratio=0.8, 
             gap.after = 1.5)

# Make circles
colourpal_nan <- colourpal[which(colourpal != "gray")]
for (i in 1:length(colourpal_nan)) {
  aname <- names(colourpal_nan[i])
  p <- ggplot() + 
    geom_circle(aes(x0=0,y0=0,r=1), fill=colourpal_nan[i], colour=NA) +
    coord_fixed() +
    theme_void()
  ggsave(p, filename = paste0("figures/paper/fig1/circle", aname,".png"),  
         bg = "transparent")
}

# Make plots of connections in brain for each node
mean_brain <- MEAN_BRAIN + t(MEAN_BRAIN)
atree <- make_null_tree()
#8, 10, 19, 21 some more didn't work
i <- 16
tree_weights <- rep(0, length(atree))
tree_weights[i] <- 1
atree <- set_tree_weights(atree, tree_weights)
A <- (make_mat(atree) > 0) * mean_brain * 
  (mean_brain >= quantile(as.vector(mean_brain), 0.0))
#A[61, 61] <- 1 # for 21 
#A[68, 68] <- 1 # for 19
#A[11, 11] <- 1 # for 8
#A[27, 27] <- 1 # for 10
#node_color <- ifelse(grepl("rh-",MAT_ROI), 2, 3)
brainconn::brainconn(atlas = "dk68", 
                     conmat = A, 
                     node.size = 3, 
                     edge.alpha = 0.6,
                     edge.color = colourpal[[NODES[i]]],
                     #node.color = node_color,
                     node.color = "gray30",
                     #view = "ortho",
                     #labels = T,
                     #all.nodes = T,
                     scale.edge.width = c(1,7),
                     show.legend = F)
ggsave(filename = paste0("figures/paper/fig1/", NODES[i],"-conn.png"),  
       bg = "transparent")

# Deskian atlas
ggplot() +
  geom_brain(atlas = dk) +
  annotate("text", x = 300, y = -50, label = "left hemisphere") +
  annotate("text", x = 1100, y = -50, label = "right hemisphere") +
  scale_fill_brain("dk") +
  theme(legend.position = "none",
        axis.line = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.background = element_blank())
ggsave(filename = paste0("figures/paper/fig1/desikan.png"),  
       bg = "transparent")
