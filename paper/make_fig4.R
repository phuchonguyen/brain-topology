# Canonical Correlation Analysis comparing performance of Tree vs. Adjacency Matrix
# representation of brain connectiomes
# Figure 4 in Tree Representations of Brain Structural Connectivity via 
# Persistent Homology by D. Li, P. Nguyen, Z. Zhang and D. Dunson
#
# Author: Phuc Nguyen, Aug. 2022

library(R.matlab)
library(rmatio)
library(tidyverse)
library(CCA)
source("tools/make_tree.R")
set.seed(123)

# Load traits
cog <- readRDS("data/CogMeasure_Aligned.Rds")
# Reverse value order for some traits in Alcohol that involves frequency
rv_cog <- c(58, 59, 60, 64, 65, 66)
for(i in rv_cog) {
  cog[,i] <- max(cog[,i], na.rm = TRUE) - cog[,i] + 1
}


# Load adjacency matrices
nets_mat <- readMat("data/matrices.mat")[[1]]
nets <- t(apply(nets_mat, 3, function(m) m[upper.tri(m)] %>% as.vector))
# Remove the columns with zero variance
id_novar <- which(apply(nets, 2, var) == 0)
nets_var <- nets[, -id_novar]

# Load regions names
ROI <- read.table("data/Desikan_Cortical_ROI_Order.txt", 
                  fill = TRUE, 
                  col.names = c("id", "name", "v1", "v2", "v3", "v4")) %>%
  filter(id > 1000) %>%
  separate(name, into = c("layer", "side", "desikan"), sep = "-") %>%
  mutate(roi = paste(side, desikan, sep = "-")) %>%
  dplyr::select(roi) %>%
  as.matrix() %>% as.vector()

# Load trees
names <- get_tree_names(make_null_tree())
trees <- readMat("data/trees.mat")[[1]]
colnames(trees) <- names
strees <- trees[, colSums(trees==0)!=1065] # remove columns with only 0 entries, 23 nodes remain

# Selecting traits
cog_pos <- c(rep("Positive", 15),
             "Negative", # Emo Rec CRTlonger time is negative
             rep("Negative", 156-150+1), # Negative emo
             rep("Positive", 160-157+1), # Positive emo
             rep("Negative", 78-74+1 + 72-68+1), # Tobacco
             rep("Negative", 88-87+1 + 84-79+1) # Illicit Drug
)
cog_id <- c(1,5,7,8,11,23,24,25,32,36,38,40, 41, 43,45,
            144, # 143:149, # Emotion Recog
            150:157, # Negative emo
            158:160, # Positive emo
            # 167:171, # Personality
            68:72, 74:78, # Tobacco
            #57:67,# Alcohol use
            79:84, 87:88 # Illicit Drug Use
)
cog_names <- c(# Cognitive skills
  "Fluid Intelligence",
  "Oral Reading Recognition",
  "Picture Vocabulary",
  "Verbal Episodic Memory",
  "Processing Speed",
  "Delay Discounting AUC 200",
  "Delay Discounting AUC 40K",
  "Spatial Orientation",
  "Sustained Attention",
  "Working Memory", # List Sort Visual
  "Picture Sequence Memory",
  "Education",
  "Income",
  "Cognitive Flexibility",
  "Visual Inhibition",  #  (Flanker)
  # Emotion Recognition
  "Emotion Recognition Time",
  # "ER40_CR","ER40ANG","ER40FEAR","ER40HAP","ER40NOE","ER40SAD",
  # Negative emotion
  "Anger-Affect", "Anger-Hostility","Anger-Physical Aggression","Fear-Affect","Fear-Somatic Arousal",
  "Sadness","Perceived Stress",
  # Positive Emotion
  "Self-Efficacy", "Life Satisfaction","Meaning-Purpose","Positive-Affect",
  # Personality
  # "Agreeableness", "Openness","Conscientiousness","Neuroticism","Extroversion",
  # Tobacco use
  "TB Age 1st_Cig","TB_DSM_Difficulty_Quitting","TB_Max_Cigs","TB_Reg_CPD","TB_Yrs_Smoked",
  "Times_Used_Any_Tobacco_Today","Avg_Weekend_Any_Tobacco_7days","Total_Cigars_7days","Avg_Weekend_Cigars_7days","Num_Days_Used_Any_Tobacco_7days",  
  # Alcohol use
  # "Alc_12_Drinks_Per_Day","Alc_12_Frq","Alc_12_Frq_5plus","Alc_12_Frq_Drk","Alc_12_Max_Drk","Alc_Age_1st_Use",
  # "Alc_Hvy_Drinks_Per_Day","Alc_Hvy_Frq","Alc_Hvy_Frq_5plus","Alc_Hvy_Frq_Drk","Alc_Hvy_Max_Drinks",
  # Illicit drug use
  "Times_Used_Illicits","Times_Used_Cocaine","Times_Used_Hallucinogens","Times_Used_Opiates","Times_Used_Sedatives",
  "Times_Used_Stimulants","Mj_Age_1st_Use","Mj_Times_Used"
)
cog_names <- gsub("_", " ", cog_names)

# CCA tree
# Canonical Correlation Analysis for Tree
z <- cog[,cog_id]
# Remove traits with more than 10% missing values
keep_id <- colMeans(is.na(z)) < 0.1 & apply(z, 2, var, na.rm=TRUE) != 0
x <- z[, keep_id] 
# Impute the missing values with means
x <- replace_na(x, as.list(colMeans(x, na.rm=T)))
colnames(x) <- cog_names[keep_id]
x <- scale(x)
# Gaussianise
zscores <- function(y){ qnorm( rank(y)/( length(y)+1 ) ) }
#x <- apply(x, 2, zscores)
y <- strees
y <- apply(y, 2, zscores)
# Canonical correlation
cct <- CCA::cc(X=x, Y=y)
cc2 <- cct$scores
# Test for significance
CCP::p.asym(cct$cor, dim(x)[1], dim(x)[2], dim(y)[2], tstat = "Wilks")
# Save result in dataframe
dat_tree <- data.frame(corr = cc2$corr.X.yscores[,1], 
                       trait = cct$names$Xnames,
                       positive = cog_pos[keep_id],
                       coef = cct$xcoef[,1]
) %>%
  mutate(sign = ifelse(corr>0, "+", "-")) %>%
  mutate(source="tree") %>%
  arrange(corr)

# Canonical Analysis for Matrix
# Demean and standardized
nets_scaled <- scale(nets_var)
# Get PCA
k <- ncol(y)
nets_svd <- svd(nets_scaled)
nets_pc <- nets_svd$u[,1:k] %*% diag(nets_svd$d[1:k])
# Gaussianise
ym <- apply(nets_pc, 2, zscores)
ccm <- CCA::cc(x, ym)
cc3 <- ccm$scores
# Correlation between two
ccm$cor[1]
# Significant test
CCP::p.asym(ccm$cor, dim(x)[1], dim(x)[2], dim(ym)[2], tstat = "Wilks")
# Save results in dataframe
dat_mat <- data.frame(corr = cc3$corr.X.yscores[,1], 
                      trait = cog_names[keep_id],
                      positive = cog_pos[keep_id],
                      coef = ccm$xcoef[,1]) %>% 
  mutate(source="mat")%>%
  arrange(corr)
# Viz
dat_tree_mat <- bind_rows(dat_tree, dat_mat) %>%
  group_by(source) %>%
  arrange(corr) %>%
  mutate(rank = rank(corr),
         rank_adj = rank + cumsum(5*abs(coef))) %>% 
  ungroup()

dat_tree_mat %>%
  ggplot() +
  labs(y="Correlation between behavioral traits \nand the first canonical variate") +
  geom_text(aes(x=source, 
                y=rank_adj,
                label=trait, color=positive,
                #size=abs(corr)
                size=coef^2)) +
  scale_x_discrete(labels= c("Adjacency Matrix Representation", "Tree Representation")) +
  scale_color_discrete(name="", labels=c("Undesirable Trait", "Desirable Trait")) +
  scale_y_continuous(breaks = c(2.9, 25, 50, 64.8),
                     labels = c(-0.60, -0.20, 0, 0.20)) +
  theme(legend.position = "bottom",
        axis.line.y = element_line(colour="gray"),
        axis.line.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.title.x = element_blank(),
        panel.background = element_blank()) +
  guides(size=FALSE)

ggsave("figures/paper/fig4/fig4-v4.eps", width = 7, height = 6.5, device=cairo_ps)
# v2 is sized by corr