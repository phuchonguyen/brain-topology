# Load pacakges
library(tidyverse)
library(BMS)
library(R.matlab)
library(viridis)
source("tools/make_tree.R")
source("tools/chordplot.R")
source("tools/tree_to_mat.R")

set.seed(1)

# Load trees
trees <- readMat("data/trees.mat")[[1]]
colnames(trees) <- NAME
strees <- trees[, colSums(trees==0)!=1065] # remove columns with only 0 entries
# Make a tree with sum and difference
reg_names <- gsub("^[lr]", "", colnames(strees))
sumdif_trees <- NULL
sum_trees <- NULL
sumdif_trees_names <- c()
sum_trees_names <- c()
sym_cor <- c()
sumdif_cor <- c()
for (n in unique(reg_names)) { # sum only lh-rh, lh-frontallobe/rh-frontallobe, lh-parietallobe/rh-parietallobe
  temp <- as.matrix(strees[, which(reg_names == n)])
  if (n %in% c("h", "h-frontallobe", "h-parietallobe")) {
    sm <- temp[,1] + temp[,-1]
    df <- temp[,1] - temp[,-1]  # left hemis - right hemis
    ad <- cbind(sm, df)
    ad2 <- sm
    sym_cor <- c(sym_cor, cor(temp[,1], temp[,-1]))
    sumdif_cor <- c(sumdif_cor, cor(sm, df))
    sumdif_trees_names <- c(sumdif_trees_names, paste(n, "lsr", sep="-") , paste(n, "lmr", sep="-"))
    sum_trees_names <- c(sum_trees_names, paste(n, "lsr", sep="-"))
  } else {
    ad <- temp
    ad2 <- temp
    if (n == "ctx") {
      sumdif_trees_names <- c(sumdif_trees_names, n)
      sum_trees_names <- c(sum_trees_names, n)
    } else {
      sumdif_trees_names <- c(sumdif_trees_names, paste0("l", n), paste0("r", n))
      sum_trees_names <- c(sum_trees_names, paste0("l", n), paste0("r", n))
    }
  }
  if (is.null(sumdif_trees)) {
    sumdif_trees <- ad
    sum_trees <- ad2
  } else {
    sumdif_trees <- cbind(sumdif_trees, ad)
    sum_trees <- cbind(sum_trees, ad2)
  }
}
colnames(sumdif_trees) <- sumdif_trees_names 
colnames(sum_trees) <- sum_trees_names
# Color map
# Given value [0,1], return an RGB value for a colour on the given ramp
map_colors <- colorRamp(viridis(50))
map_colors <- colorRamp(c("#0072B2",
                          "#D3D3D3", #gray #"#FFFFFF", # white
                          "#D55E00"))

# Scale values
my_scale <- function(tree, type="unit", tree.max=NULL, tree.min=NULL) {
  if (is.null(tree.max)) {tree.max = max(tree)}
  if (is.null(tree.min)) {tree.min = min(tree)}
  if (type=="unit") {
    if (tree.max - tree.min == 0) return(tree - tree.min)
    return((tree - tree.min)/(tree.max - tree.min))
  }
  if (type=="polar") {
    if (max(abs(tree.min), abs(tree.max)) == 0) return(tree + 0.5)
    return(tree/(2*max(abs(tree.min), abs(tree.max))) + 0.5)
  }
}

# Load traits
cog <- readRDS("data/CogMeasure_Aligned.Rds")

# Helper function
# normalize
zscores <- function(y){ qnorm( rank(y)/( length(y)+1 ) ) }
# format coef into tree
make_coef_tree <- function(coef_select) {
  nulltree <- make_null_tree()
  for (i in 1:length(coef_select)) {
    n <- names(coef_select)[i]
    if (n == "(Intercept)") next
    if (!grepl("h-lsr|h-lmr|h-frontallobe|h-parietallobe", n)) {
      nulltree[n] <- coef_select[n]
    } else {
      n <- gsub("\\.", "-", n)
      if (str_detect(n, "lsr")) {
        n <- gsub("-lsr", "", n)
        n <- paste0(c("r", "l"),n)
        nulltree[which(names %in% n)] <- coef_select[i]
      } else if (str_detect(n, "lmr")) {
        neg <- 1*(coef_select[i] < 0)
        n <- gsub("-lmr", "", n)
        if (neg) {
          n <- paste0(c("r"),n)
          nulltree[n] <- nulltree[n] -1*coef_select[i]
        } else {
          n <- paste0(c("l"),n)
          nulltree[n] <- nulltree[n] + coef_select[i]
        }
      } else print(n)
    }
  }
  return(nulltree)
}


fit_bms_yid <- function(y_id, use_sum=TRUE, scale=F, nmodel=2000){
  y <- cog[,y_id] #%>% zscores()
  if (use_sum) {bms_dat <- cbind(y, sumdif_trees)}
  else {bms_dat <- cbind(y, strees)}
  if (scale) bms_dat <- scale(bms_dat)
  bms_dat <- bms_dat[complete.cases(bms_dat),]
  mod <- bms(bms_dat, 
             burn = 50000, iter = 1e+05, g = "UIP",
             mprior = "uniform", nmodel = nmodel, mcmc = "bd", user.int = F)
  return(mod)
}

plot_bms <- function(mod, title="", tree.min=NULL, tree.max=NULL,
                     consider.se=FALSE) {
  coef_select <- coef(mod, condi.coef=FALSE)[,2]#[coef(list_mod)[,1] > 0.5] # inclusion prob > 50%
  coef_pip <- coef(mod)[,1]
  if (consider.se) {
    not_sig <- (sign(coef(mod)[,2]) != sign(coef(mod)[,2]+1.96*coef(mod)[,3])) | (sign(coef(mod)[,2]) != sign(coef(mod)[,2]-1.96*coef(mod)[,3]))
    coef_select[not_sig] <- 0
  }
  A_tree <- make_coef_tree(coef_select)
  A_tree_scaled <- my_scale(A_tree, type="polar", tree.min=tree.min, tree.max=tree.max)
  A_tree_pip <- make_coef_tree(coef_pip)
  colorpal <- map_colors(A_tree_scaled)
  colorpal <- apply(cbind(colorpal, A_tree_pip), 1, function(x) rgb(x[1], x[2], x[3],
                                                                    alpha=min(255, 225*x[4] + 30),
                                                                    maxColorValue=255))
  # colorpal[A_tree == 0] <- "blank"
  # colorpal[A_tree == 0] <- rgb(211,211,211,alpha=50,maxColorValue=255)
  names(colorpal) <- names(A_tree)
  chord_tree(ribcol = colorpal, legend = T, tree = A_tree, 
             title = title,
             palette = "redwhiteblue", legend.x = 1.05, 
             img.labeled = T)
}

plot_bms_nosum <- function(mod, title="", tree.min=NULL, tree.max=NULL, 
                           pip.thres=1, consider.se=F) {
  coef_estimate <- coef(mod, condi.coef=FALSE)[,2]
  if (consider.se) {
    not_sig <- (sign(coef(mod)[,2]) != sign(coef(mod)[,2]+1.96*coef(mod)[,3])) | (sign(coef(mod)[,2]) != sign(coef(mod)[,2]-1.96*coef(mod)[,3]))
    coef_estimate[not_sig] <- 0
  }
  coef_estimate[coef(mod)[,5]] <- coef_estimate # reorder
  coef_pip <- coef(mod)[,1]
  coef_pip[coef(mod)[,5]] <- coef_pip
  # NEW
  coef_estimate <- coef_estimate * (coef_pip >= pip.thres)
  names(coef_estimate) <- names(coef_pip) <- colnames(strees)
  A_tree_scaled <- my_scale(coef_estimate, type="polar", tree.min=tree.min, tree.max=tree.max)
  A_tree_pip <- coef_pip #ifelse(coef_pip >= pip.thres, 1, coef_pip)
  colorpal <- map_colors(A_tree_scaled)
  colorpal <- apply(cbind(colorpal, A_tree_pip), 1, function(x) rgb(x[1], x[2], x[3],
                                                                    alpha=min(255, 225*x[4] + 30),
                                                                    maxColorValue=255))
  # colorpal[A_tree == 0] <- "blank"
  # colorpal[A_tree == 0] <- rgb(211,211,211,alpha=50,maxColorValue=255)
  names(colorpal) <- names(A_tree_scaled)
  chord_tree(ribcol = colorpal, legend = T, tree = coef_estimate, 
             title = title,
             palette = "redwhiteblue", legend.x = 1.05, 
             img.labeled = T)
}

transparent_bg <- rgb(1, 1, 1, alpha=0)
PIP_THRES <- 0.75
CONSIDER_SE <- F
# Working Memory (Visual / List Sort)
# This task assesses working memory and requires the participant to sequence 
# different visually- and orally-presented stimuli. Pictures of different foods 
# and animals are displayed with both a sound clip and written text that name the item. 
# The task has two different conditions: 1-List and 2-List. 
# In the 1-List condition, participants are required to order a series of objects 
# (either food or animals) in size order from smallest to largest. 
# In the 2-List condition, participants are presented both food and animals
# and are asked to report the food in size order, followed by the animals in size order. 
# Children ages 3-6 have four practice items in each condition: two practice items 
# in which the images appear simultaneously on the screen and two practice items in 
# which the images briefly "flash" sequentially on the screen. 
# Participants ages 7-85 have two practice items, both "flashing" in each condition. 
# Different instructions are provided for 3-6 and 7-85 year olds in English and for 3-6, 7-17 and 18-85 year olds in Spanish.
y_id <- 36
# tree
wm_mod <- fit_bms_yid(36, use_sum = F, scale = T)
svglite::svglite(paste0("figures/fig7/Working Memory ", PIP_THRES, CONSIDER_SE, ".svg"), 
                 width = 7, height = 5.5,
                 bg=transparent_bg)
plot_bms_nosum(wm_mod, title = "Working Memory", 
               consider.se = CONSIDER_SE, pip.thres = PIP_THRES)
dev.off()

# Visual Inhibitor
# The Flanker task measures both a participant's attention and inhibitory control. 
# The test requires the participant to focus on a given stimulus while inhibiting 
# attention to stimuli (fish for ages 3-7 or arrows for ages 8-85) flanking it. 
# Sometimes the middle stimulus is pointing in the same direction as the "flankers"
# (congruent) and sometimes in the opposite direction (incongruent). 
# Scoring is based on a combination of accuracy and reaction time, 
# and the test takes approximately 3 minutes to administer. 
# This test is recommended for ages 3-85.
y_id <- 45
vi_mod <- fit_bms_yid(45, use_sum = F, scale = T)
svglite::svglite(paste0("figures/fig7/Visual Inhibitor ", PIP_THRES, CONSIDER_SE, ".svg"), 
                 width = 7, height = 5.5,
                 bg=transparent_bg)
#plot_bms(fit_bms_yid(y_id), title = "Visual Inhibitor")
plot_bms_nosum(vi_mod, title = "Visual Inhibitor", 
               consider.se = CONSIDER_SE, pip.thres = PIP_THRES)
dev.off()


# Fluid Intel
# Fluid intelligence is measured using Raven’s Progressive Matrices 
# (Prabhakaran et al. 1997; Christoff et al. 2001; Gray et al. 2003; Conway et al. 2005;
# Gray et al. 2005; Wendelken et al. 2008). 
# We use Form A of an abbreviated version of the Raven’s developed by Gur
# and colleagues (Bilker et al. 2012). 
# Participants are presented with patterns made up of 2x2, 3x3 or 1x5 arrangements of squares,
# with one of the squares missing. The participant must pick one of five response
# choices that best fits the missing square on the pattern. 
# The task has 24 items and 3 bonus items, arranged in order of increasing difficulty. 
# However, the task discontinues if the participant makes 5 incorrect responses in a row.
y_id <- 1
fi_mod <- fit_bms_yid(1, use_sum = F, scale = T)
svglite::svglite(paste0("figures/fig7/Tree Fluid Intelligence ", PIP_THRES, CONSIDER_SE, ".svg"),
                 width = 7, height = 5.5,
                 bg=transparent_bg)
plot_bms_nosum(fi_mod, title = "Fluid Intelligence", 
               consider.se = CONSIDER_SE, pip.thres = PIP_THRES)
# plot_bms_nosum(fi_mod, title = "Fluid Intelligence", pip.thres = 1.0,
#                consider.se = T)
dev.off()


# Spatial Orientation
# Spatial orientation processing is measured using the Variable Short Penn Line Orientation
# Test (Gur et al. 2001a; Gur et al. 2010). 
# Participants are shown two lines with different orientations. 
# They have to rotate one of the lines (a moveable blue one) so that is 
# parallel to the other line (a fixed red line). 
# The rotation of the blue line is accomplished by clicking buttons on 
# the keyboard that rotate the lines either clockwise or counterclockwise. 
# Across trials, the lines vary in their relative location on the screen, 
# though the distance between the centers of the two lines is always the same. 
# The length of the red line is always the same, but the length of the blue line 
# can be either short or long. 
# There are a total of 24 trials
y_id <- 25
so_mod <- fit_bms_yid(25, use_sum = F, scale = T, nmodel = 10000)
svglite::svglite(paste0("figures/fig7/Tree Spatial Orientation ", PIP_THRES, CONSIDER_SE, ".svg"),
                 width = 7, height = 5.5,
                 bg=transparent_bg)
plot_bms_nosum(so_mod, title = "Spatial Orientation",
               consider.se = CONSIDER_SE, pip.thres = PIP_THRES)
# plot_bms_nosum(so_mod, title = "Spatial Orientation", pip.thres = 1,
#                consider.se = T)
dev.off()


# Oral Reading
# Separate but parallel reading tests have been developed in English and in Spanish. 
# In either language, the participant is asked to read and pronounce letters 
# and words as accurately as possible. The test administrator scores them as right or wrong. 
# For the youngest children, the initial items require them to identify letters 
# (as opposed to symbols) and to identify a specific letter in an array of 4 symbols. 
# The test is given via a computerized adaptive format and requires approximately 3 minutes. 
# This test is recommended for ages 7-85, but is available for use as young as age 3, if requested.
y_id <- 5
or_mod <- fit_bms_yid(5, use_sum = F, scale = T)
svglite::svglite(paste0("figures/fig7/Oral Reading ", PIP_THRES, CONSIDER_SE, ".svg"), 
                 width = 7, height = 5.5,
                 bg=transparent_bg)
plot_bms_nosum(or_mod, title = "Oral Reading Recognition", 
               consider.se = CONSIDER_SE, pip.thres = PIP_THRES)
# plot_bms_nosum(or_mod, title = "Oral Reading Recognition", pip.thres = 1,
#                consider.se = T)
dev.off()


# Pic Vocab
# This measure of receptive vocabulary is administered in a computerized adaptive format.
# The respondent is presented with an audio recording of a word and four photographic 
# images on the computer screen and is asked to select the picture that most closely 
# matches the meaning of the word. 
# This test takes approximately 4 minutes to administer and is recommended for ages 3-85.
y_id <- 7
pv_mod <- fit_bms_yid(7, use_sum = F, scale = T)
svglite::svglite(paste0("figures/fig7/Picture Vocabulary ", PIP_THRES, CONSIDER_SE, ".svg"),
                 width = 7, height = 5.5,
                 bg=transparent_bg)
#plot_bms(fit_bms_yid(y_id), title = "Picture Vocabulary")
plot_bms_nosum(pv_mod, title = "Picture Vocabulary", 
               consider.se = CONSIDER_SE, pip.thres = PIP_THRES)
dev.off()



