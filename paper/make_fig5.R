# Make figure to compare predictive performance of Tree vs. Adjacency Matrix
# representation of brain connectomes.
# Figure 5 and 6 in Tree Representations of Brain Structural Connectivity via 
# Persistent Homology by D. Li, P. Nguyen, Z. Zhang and D. Dunson
#
# Author: Phuc Nguyen, Aug. 2022

library(ggplot2)
library(ggrepel)
library(tidyverse)
library(gridExtra)

abbre_df <- data.frame(acronym=c("FI", "ORR", "PV", "SO"),
                       trait=c("Fluid\nintelligence", "Oral reading\nrecognition", "Picture\nvocabulary", "Spatial\norientation")) 
cog_names <- c(# Cognitive skills
  paste("FI", c("correct responses", "skipped items", "median reaction time")),
  paste("ORR", c("unadjusted", "age-adjusted")),
  paste("PV", c("unadjusted", "age-adjusted")),
  rep("Verbal Episodic Memory", 2),
  rep("Processing Speed", 2),
  rep("Delay Discounting AUC", 13),
  #"Delay Discounting AUC 40K",
  paste("SO", c("number correct", "median reaction time", "number off")),
  rep("Sustained Attention", 8),
  rep("Working Memory", 2), # List Sort Visual
  rep("Picture Sequence Memory", 2),
  "Education",
  "Income",
  rep("Cognitive Flexibility", 2),
  rep("Visual Inhibition", 2)  #  (Flanker)
)
cog_abbre <- c(
  rep("FI - Fluid intelligence", 3),
  rep("ORR - Oral reading recognition", 2),
  rep("PV - Picture vocabulary", 2),
  rep("Verbal Episodic Memory", 2),
  rep("Processing Speed", 2),
  rep("Delay Discounting AUC", 13),
  #"Delay Discounting AUC 40K",
  rep("SO - Spatial orientation", 3),
  rep("Sustained Attention", 8),
  rep("Working Memory", 2), # List Sort Visual
  rep("Picture Sequence Memory", 2),
  "Education",
  "Income",
  rep("Cognitive Flexibility", 2),
  rep("Visual Inhibition", 2)  #  (Flanker)
)
cog_id <- c(1,5,7,8,11,23,24,25,32,36,38,40, 41, 43,45)
names <- rep("", 45)
names[1:7] <- cog_names[1:7]
data_5a <- read.csv("data/Fig5_1.csv", header = FALSE) %>%
  t() %>%
  as.data.frame() %>%
  magrittr::set_colnames(c("matrix", "tree")) %>%
  mutate(name = cog_names)
p5a <- ggplot(data_5a, aes(x=100*matrix, y=100*tree, label=name)) + 
  geom_point() +
  geom_point(data = subset(data_5a, tree > 0.025 | matrix > 0.025),
             aes(x=100*matrix, y=100*tree), colour="brown3", size=1.5) +
  geom_text_repel(
    data = subset(data_5a, tree > 0.025 & tree > matrix),
    nudge_y       = c(-0.5, 0.75, 1.5, 1),
    nudge_x       = -c(5, 4, 3, 3.5),
    box.padding   = 0.5,
    point.padding = 0.5,
    force         = 0,
    segment.size  = 0.2,
    segment.color = "grey75",
    max.overlaps  = 50,
    direction     = "both"
    ) +
  geom_text_repel(
    data = subset(data_5a, matrix > 0.05 & tree < matrix),
    nudge_y       = - c(-2.5, 5, 2),
    nudge_x       = - c(2.5, 0, 0.25),
    box.padding   = 0.5,
    point.padding = 0.5,
    force         = 0,
    segment.size  = 0.2,
    segment.color = "grey75",
    max.overlaps = 50,
    direction     = "both"
  ) +
  geom_abline(slope=1, intercept = 0) +
  xlim(-10, 10) + ylim(-10, 10) +
  theme_light(base_size = 18) +
  # theme(axis.title=element_text(size=14)) +
  labs(x="Adjacency matrix (% change)", y="Tree (% change)", title="Linear Regression")
ggsave(plot = p5a, filename = "figures/fig5/fig5a.pdf", device="pdf", height = 5, width = 5)

data_5b <- read.csv("data/Fig5_2.csv", header = FALSE) %>%
  t() %>%
  as.data.frame() %>%
  magrittr::set_colnames(c("matrix", "tree")) %>%
  mutate(name = cog_names,
         abbre = cog_abbre)
p5b <- ggplot(data_5b, aes(x=100*matrix, y=100*tree, label=name)) + 
  geom_point() + 
  geom_point(data = subset(data_5b, tree > 0.05 | matrix > 0.05),
             aes(x=100*matrix, y=100*tree), colour="brown3", size=1.5) +
  geom_text_repel( 
    data = subset(data_5b, tree > 0.05 & tree > matrix),
    nudge_y       = c(-1., 0.75, 0.75, 0.5, 1),
    nudge_x       = - c(2, 3.5, 2, 3.5, 1),
    box.padding   = 0.5,
    point.padding = 0.5,
    force         = 0,
    segment.size  = 0.2,
    segment.color = "grey75",
    max.overlaps = 50,
    direction     = "both"
  ) +
  geom_text_repel( 
    data = subset(data_5b, matrix > 0.05 & tree < matrix),
    nudge_y       = - c(-1.5, 2.5, 1.5),
    nudge_x       = - c(1.5, -0.5, 0.25),
    box.padding   = 0.5,
    point.padding = 0.5,
    force         = 0,
    segment.size  = 0.2,
    segment.color = "grey75",
    max.overlaps = 50,
    direction     = "both"
  ) +
  geom_abline(slope=1, intercept = 0) +
  xlim(-5, 10) + ylim(-5, 10) +
  theme_light(base_size = 18) +
  labs(x="Adjacency matrix (% change)", y="Tree (% change)", title="Gaussian Process Regression")
ggsave(plot = p5b, filename = "figures/fig5/fig5b.pdf", device="pdf", height = 5, width = 5)


p5 <- grid.arrange(p5a, p5b, 
                   tableGrob(abbre_df, rows = NULL, theme=ttheme_default()), 
                   nrow = 1, widths = c(5, 5, 2))
ggsave(plot = p5, filename = "figures/fig5/fig5.pdf", device="pdf", height = 5, width = 12)

data_6a <- read.csv("data/Fig6_1.csv", header = FALSE) %>%
  t() %>%
  as.data.frame() %>%
  magrittr::set_colnames(c("matrix", "tree")) %>%
  mutate(name = cog_names)
p6a <- ggplot(data_6a, aes(x=matrix, y=tree, label=name)) + 
  geom_point() + 
  geom_text_repel( 
    data = subset(data_6a, tree > 0.2 & tree > matrix),
    nudge_y       = -c(0.0, -0.01, -0.005, -0.02),
    nudge_x       = -c(.08, .06, .08, .03),
    box.padding   = 0.5,
    point.padding = 0.5,
    force         = 0,
    segment.size  = 0.2,
    segment.color = "grey75",
    max.overlaps = 50,
    direction     = "both"
  ) +
  geom_text_repel( 
    data = subset(data_6a, matrix > 0.2 & tree < matrix),
    nudge_y       = - c(0.06, 0.08, 0.11),
    nudge_x       = - c(0., 0.0, .0),
    box.padding   = 0.5,
    point.padding = 0.5,
    force         = 0,
    segment.size  = 0.2,
    segment.color = "grey75",
    max.overlaps = 50,
    direction     = "both"
  ) +
  geom_abline(slope=1, intercept = 0) +
  xlim(-0.05, 0.3) + ylim(-0.05, 0.3) +
  theme_light(base_size = 18) +
  labs(x="Adjacency matrix (correlation)", y="Tree (correlation)", title="Linear Regression")
ggsave(plot=p6a, filename = "figures/fig5/fig6a.pdf", device="pdf", 
       height = 5, width = 5)


data_6b <- read.csv("data/Fig6_2.csv", header = FALSE) %>%
  t() %>%
  as.data.frame() %>%
  magrittr::set_colnames(c("matrix", "tree")) %>%
  mutate(name = cog_names)
p6b <- ggplot(data_6b, aes(x=matrix, y=tree, label=name)) + 
  geom_point() + 
  geom_text_repel( 
    data = subset(data_6b, tree > 0.2 & tree > matrix),
    nudge_y       = -c(0.0, 0.015, -0.015, -0.01, -0.02),
    nudge_x       = -c(.1, .05, .07, .07, .05),
    box.padding   = 0.5,
    point.padding = 0.5,
    force         = 0,
    segment.size  = 0.2,
    segment.color = "grey75",
    max.overlaps = 50,
    direction     = "both"
  ) +
  geom_text_repel( 
    data = subset(data_6b, matrix > 0.2 & tree < matrix),
    nudge_y       = - c(-0.05, 0.03, 0.05),
    nudge_x       = - c(0.03, 0.0, .0),
    box.padding   = 0.5,
    point.padding = 0.5,
    force         = 0,
    segment.size  = 0.2,
    segment.color = "grey75",
    max.overlaps = 50,
    direction     = "both"
  ) +
  geom_abline(slope=1, intercept = 0) +
  xlim(-0.05, 0.35) + ylim(-0.05, 0.35) +
  theme_light(base_size = 18) +
  labs(x="Adjacency matrix (correlation)", y="Tree (correlation)", title="Gaussian Process Regression")
ggsave(plot=p6b, filename = "figures/fig5/fig6b.pdf", device="pdf", 
       height = 5, width = 5)

p6 <- grid.arrange(p6a, p6b, 
                   tableGrob(abbre_df, rows = NULL, theme=ttheme_default()), 
                   nrow = 1, widths = c(5, 5, 2))
ggsave(plot = p6, filename = "figures/fig5/fig6.pdf", device="pdf", height = 5, width = 12)
