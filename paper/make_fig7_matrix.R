# Regress traits on PCA of vectorized adjacency matrix
library(R.matlab)
library(rmatio)
library(tidyverse)
library(elasticnet)
library(glmnet)
library(caret)
library(CCA)
library(BMS)
source("tools/make_tree.R")
source("tools/chordplot.R")
source("tools/tree_to_mat.R")

set.seed(123)

# Load adjacency matrices
nets_mat <- read.mat("data/HCP_cortical_TensorData_desikan.mat")[["loaded_tensor"]][,,1,]
nets <- t(apply(nets_mat, 3, as.vector))
# Remove the columns with zero variance
id_novar <- which(apply(nets, 2, var) == 0)
nets_var <- nets[, -id_novar]
# Demean and standardize
nets_scaled <- scale(nets_var)

# Load regions names
ROI <- read.table("data/Desikan_Cortical_ROI_Order.txt", 
                  fill = TRUE, 
                  col.names = c("id", "name", "v1", "v2", "v3", "v4")) %>%
  filter(id > 1000) %>%
  separate(name, into = c("layer", "side", "desikan"), sep = "-") %>%
  mutate(roi = paste(side, desikan, sep = "-")) %>%
  dplyr::select(roi) %>%
  as.matrix() %>% as.vector()
ROI_var <- ROI[-id_novar]

# Load traits
cog <- readRDS("data/CogMeasure_Aligned.Rds")

cog_id <- c(1,5,7,8,11,23,24,25,34,36,38,43,45
)
cog_names <- c("Fluid Intelligence",
               "Oral Reading Recognition",
               "Picture Vocabulary",
               "Verbal Episodic Memory",
               "Processing Speed",
               "Delay Discounting AUC 200",
               "Delay Discounting AUC 40K",
               "Spatial Orientation",
               "Sustained Attention",
               "Visual Memory (list sort)",
               "Picture Sequence Memory",
               "Cognitive Flexibility",
               "Visual Inhibition (Flanker)"
)
cog_id <- c(1,5,7, 25,36,45)
cog_names <- c("Fluid Intelligence",
               "Oral Reading Recognition",
               "Picture Vocabulary",
               "Spatial Orientation",
               "Working Memory",
               "Visual Inhibitor"
               )

# Source: http://www2.stat.duke.edu/~pdh10/Teaching/832/Notes/copula.pdf
zscores <-function(y,ties.method="average")
{
  u<-rank(y,na.last="keep",ties.method=ties.method)/(sum(!is.na(y))+1)
  z<-qnorm(u)
  names(z)<-names(y)
  m<-dim(y)
  if(length(m)==2){
    z<-matrix(z,nrow=m[1],ncol=m[2])
    dimnames(z)<-dimnames(y)
  }
  if(length(m)>=3){
    z<-array(z,dim=m)
    dimnames(z)<-dimnames(y)
  }
  z
}
# Do PCA
CONSIDER_SE <- F
ALPHA <- 0.05 #/ length(cog_id)  # correct for multiple testing NOT PLOTTED
PIP_THRES <- 0.0
K_ops <- c(23)
n_ops <- c(50)
udv <- svd(nets_scaled)
transparent_bg <- rgb(1, 1, 1, alpha=0)
for (j in 1:length(cog_id)) {
  i <- cog_id[j]
  iname <- cog_names[j]
  
  for (k in K_ops) {
    pcs <- udv$u[,1:k] %*% diag(udv$d[1:k])
    # temp_df <- cbind(scale(cog[,i]), pcs) %>% 
    #   as.data.frame() %>%
    #   magrittr::set_colnames(c("y", paste0("pc", 1:23))) %>%
    #   filter(complete.cases(.))
    # Fit linear model
    #----------------------------
    # lm_mod <- lm(y ~ . - 1, temp_df)
    # ridge_cv <- cv.glmnet(x = as.matrix(temp_df[,-1]), 
    #                       y = temp_df$y, family = "gaussian", 
    #                       alpha = 0, intercept = F)
    # ridge_mod <- glmnet(x = as.matrix(temp_df[,-1]), 
    #                     y = temp_df$y, family = "gaussian", 
    #                     alpha = 0, intercept = F, lambda = ridge_cv$lambda.min)
    # #Get the coefficients and p-values without the intercept
    # lm_coef <- coef(lm_mod)
    # lm_pval <- summary(lm_mod)$coef[,4]
    # lm_sigcoef <- lm_coef
    # lm_sigcoef[lm_pval > ALPHA] <- 0
    # ridge_coef <- coef(ridge_mod)[-1]
    # Fit BMS model
    #----------------------------
    bms_dat <- cbind(cog[,i], pcs)
    bms_dat <- bms_dat[complete.cases(bms_dat), ]
    bms_mod <- bms(bms_dat,
    burn = 50000, iter = 1e+05, g = "UIP",
    mprior = "uniform", nmodel = 2000, mcmc = "bd", user.int = F)
    #Get the coefficients with >0.5 inclusion probability
    bms_sigcoef  <- coef(bms_mod)[,2]
    # if (CONSIDER_SE) {
    #   not_sig <- (sign(coef(bms_mod)[,2]) != sign(coef(bms_mod)[,2]+1.96*coef(bms_mod)[,3])) | (sign(coef(bms_mod)[,2]) != sign(coef(mod)[,2]-1.96*coef(bms_mod)[,3]))
    #   bms_sigcoef[not_sig] <- 0
    # }
    bms_sigcoef[coef(bms_mod)[,1] < PIP_THRES] <- 0 # zero out features with low PIP
    bms_sigcoef[coef(bms_mod)[,5]] <- bms_sigcoef # reorder coefs to match input
    #----------------------------
    # Map PCs coefficients back to original coefficients
    temp <- udv$v[,1:k] %*% bms_sigcoef
    nets_coef <- rep(0, ncol(nets))
    nets_coef[-id_novar] <- temp
    amat <- matrix(nets_coef, 68, 68)
    colnames(amat) <- ROI
    rownames(amat) <- ROI
    # topn <- sum(temp>2*sd(temp)+mean(temp) | temp<mean(temp)-2*sd(temp))
    for (n in n_ops) {
      svglite::svglite(paste("figures/fig7/PCA ", iname, 
                              PIP_THRES, CONSIDER_SE, n, ".svg"), 
                       width = 7.5, height = 5.9,
                       bg=transparent_bg)
      chord_adjmat(A=amat, topn=n, palette="redwhiteblue",
                   img.labeled=TRUE, legend.x = 1.05,
                   title=iname, legend=TRUE, lwd=0.5)
      dev.off()
    }
  } 
}


# ### Each model individually ###
# pcs <- udv$u[,1:k] %*% diag(udv$d[1:k])
# 
# fit_bms_pca <- function(i) {
#   bms_dat <- cbind(cog[,i], pcs)
#   bms_mod <- bms(bms_dat,
#                  burn = 50000, iter = 1e+05, g = "UIP",
#                  mprior = "uniform", nmodel = 2000, mcmc = "bd", user.int = F)
#   return(bms_mod)
# }
# 
# plot_bms_pca <- function(bms_mod, topn=50, pip.thres=0.5, title="") {
#   #Get the coefficients with >0.5 inclusion probability
#   bms_sigcoef <- coef(bms_mod)[,2]
#   bms_sigcoef[coef(bms_mod)[,1] < pip.thres] <- 0
#   bms_sigcoef[coef(bms_mod)[,5]] <- bms_sigcoef
#   # Map PCs coefficients back to original coefficients
#   temp <- udv$v[,1:k] %*% bms_sigcoef
#   nets_coef <- rep(0, ncol(nets))
#   nets_coef[-id_novar] <- temp
#   amat <- matrix(nets_coef, 68, 68)
#   colnames(amat) <- ROI
#   rownames(amat) <- ROI
#   chord_adjmat(A=amat, topn=topn, palette="redwhiteblue",
#                img.labeled=TRUE, legend.x = 1.05,
#                title=title, legend=TRUE)
# }
# 
# #Fluid Intel
# fi_pca <- fit_bms_pca(1)
# plot_bms_pca(fi_pca, topn = 100, pip.thres = 0.85, title = "Fluid Intelligence")
# 
# #Spatial Orientaiton
# so_pca <- fit_bms_pca(25)
# plot_bms_pca(so_pca, topn = 100, pip.thres = 0.85, title = "Spatial Orientation")
# 
# #Picture Vocab
# pv_pca <- fit_bms_pca(7)
# plot_bms_pca(pv_pca, topn = 100, pip.thres = 0.85, title = "Picture Vocabulary")
# 
# #Oral Reading
# or_pca <- fit_bms_pca(5)
# plot_bms_pca(or_pca, topn = 50, pip.thres = 0.9, title = "Oral Reading")
# 
# #Working Memory
# wm_pca <- fit_bms_pca(36)
# plot_bms_pca(wm_pca, topn = 100, pip.thres = 0.85, title = "Working Memory")
# 




