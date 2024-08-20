# This script creates Fig 2, Fig S6
# Please cite "Recurrent evolution of small body size and loss of the sword ornament in Northern swordtail fish" Preising et al. 2024  Evolution if you use this script

library(tidyverse)
library(ape)
library(phytools)
library(caper)
library(ggrepel)
library(ggridges)
library(diptest)
library(cowplot)
library(ggtree)
library(ggtreeExtra)
library(treeio)
library(ggbiplot)

# https://github.com/YuLab-SMU/enrichplot/issues/249
# for issues with clade labels
nodeid.tbl_tree <- utils::getFromNamespace("nodeid.tbl_tree", "tidytree")
rootnode.tbl_tree <- utils::getFromNamespace("rootnode.tbl_tree", "tidytree")
offspring.tbl_tree <- utils::getFromNamespace("offspring.tbl_tree", "tidytree")
offspring.tbl_tree_item <- utils::getFromNamespace(".offspring.tbl_tree_item", "tidytree")
child.tbl_tree <- utils::getFromNamespace("child.tbl_tree", "tidytree")
parent.tbl_tree <- utils::getFromNamespace("parent.tbl_tree", "tidytree")

# read in data
xipho_phen <- read.csv("data/phy_test_info_20240430_recode.csv")
# make vector of allometric traits to iterate over
res_traits <- c("swd_length","upper_sword_edge_width","upper_sword_edge_length","lower_sword_edge_width","lower_sword_edge_length",
                "dorsal_fin_length","dorsal_fin_height","body_depth","peduncle_depth","caudal_fin_height","caudal_fin_length",
                "peduncle_edge","lower_edge_2_peduncle","vertical_bars")

# compute residuals for allometric traits and store as vars in initial df
for (trait in res_traits){
  
  formula <- as.formula(paste(trait, "~ std_length"))
  lm_model <- lm(formula,data=xipho_phen)
  residuals <- lm_model$residuals
  varname <- paste(trait, "_residuals", sep = "")
  
  xipho_phen <- xipho_phen %>%
    dplyr::mutate(!!varname := residuals)
  
}

# store a copy of individual values for non-phylogenetic PCA later
indiv_m_phen <- xipho_phen

# figure out cutoff for sword_presence variable - max Xbir swd_index
xipho_phen %>%
  filter(species == "Xbirchmanni") %>%
  mutate(swd_index = swd_length/std_length) %>%
  group_by(species,sex) %>%
  summarise(max_swd_length = max(swd_index, na.rm=TRUE))

xipho_phen <- xipho_phen %>%
  # remove monticolous and pjonesii
  filter(!species %in% c("Xmonticolus", "Pjonesii")) %>%
  # recode sex as 0/1
  mutate(sex = case_when(sex == "M" ~ 0,
                         sex == "F" ~ 1),
         # compute sword presence variable
         sword_index = swd_length/std_length,
         sword_presence = case_when(sword_index > 0.01943906 ~ 1,
                                    sword_index <= 0.01943906 ~ 0)) %>%
  # group by grouping vars
  group_by(species,sex,tactic,clade,tree_lab) %>%
  summarise_at(vars(std_length,
                    swd_length,
                    sword_presence,
                    sword_pigmentation,
                    dorsal_fin_pigmentation,
                    vertical_bars,
                    caudal_fin_pigmentation,
                    dorsal_fin_pigmentation,
                    upper_sword_edge_width,
                    upper_sword_edge_length,
                    lower_sword_edge_width,
                    lower_sword_edge_length,
                    dorsal_fin_length,
                    dorsal_fin_height,
                    body_depth,
                    peduncle_depth,
                    caudal_fin_length,
                    caudal_fin_height,
                    false_gravid_spot,
                    peduncle_edge,
                    lower_edge_2_peduncle,
                    sc,
                    atro,
                    carbo,
                    cb,
                    false_gravid_spot,
                    swd_length_residuals,
                    upper_sword_edge_width_residuals,
                    upper_sword_edge_length_residuals,
                    lower_sword_edge_width_residuals,
                    lower_sword_edge_length_residuals,
                    lower_sword_edge_length_residuals,
                    dorsal_fin_length_residuals,
                    dorsal_fin_height_residuals,
                    body_depth_residuals,
                    peduncle_depth_residuals,
                    caudal_fin_height_residuals,
                    caudal_fin_length_residuals,
                    peduncle_edge_residuals,
                    lower_edge_2_peduncle_residuals,
                    vertical_bars_residuals),mean) %>%
  mutate(tactic = case_when(tactic == "sneaker" ~ "coercive",
                            tactic == "courter" ~ "courtship",
                            tactic %in% c("ambiguous","polymorphic") ~ tactic))

#####################
# read in tree data #
#####################

# read in tree
xipho_tree <- read.nexus("data/RAxML_bipartitions.NS_complete_20231215_treefile_addmorphs.txt.tre")
# remove kallmani
xipho_tree <- drop.tip(xipho_tree, c("Xkal"))

#################
# male analysis #
#################
# create male-only df
xipho_phen_m <- xipho_phen %>%
  filter(sex == 0) %>%
  ungroup() %>%
  dplyr::select(-c(sex,
                   swd_length,
                   vertical_bars,
                   caudal_fin_pigmentation,
                   dorsal_fin_pigmentation,
                   upper_sword_edge_width,
                   upper_sword_edge_length,
                   lower_sword_edge_width,
                   lower_sword_edge_length,
                   dorsal_fin_length,
                   dorsal_fin_height,
                   body_depth,
                   peduncle_depth,
                   caudal_fin_length,
                   caudal_fin_height,
                   peduncle_edge,
                   lower_edge_2_peduncle))
rownames(xipho_phen_m) <- xipho_phen_m$tree_lab

# check that rownames are all in the tree tip labels
rownames(xipho_phen_m) %in% xipho_tree$tip.label

# make row names in phenotype data match order of tree tips
xipho_phen_m <- xipho_phen_m %>%
  arrange(factor(tree_lab, levels = xipho_tree$tip.label)) %>%
  dplyr::select(-c(sc)) %>% #remove spotted caudal (absent in males)
  mutate(label = case_when(tree_lab %in% c("Xnig_small","Xnig_large","Xmul_small","Xmul_large","Xcon","Xmon","Xpyg","Xbir","Xhel") ~ tree_lab,
                           !(tree_lab %in% c("Xnig_small","Xnig_large","Xmul_small","Xmul_large","Xcon","Xmon","Xpyg","Xbir","Xhel")) ~ "")) %>%
  mutate(label=str_replace(label,"_.*",""))
rownames(xipho_phen_m) <- xipho_phen_m$tree_lab

#write.csv(xipho_phen_m, file="data/xipho_phen_m.csv", quote = F, row.names = F)

# keep one instance of this data as a matrix with rownames for things in phytools
xipho_phen_m_mat <- as.matrix(xipho_phen_m[5:25])
rownames(xipho_phen_m_mat) <- rownames(xipho_phen_m)

# run phylogenetic PCA
xipho_phen_m_pPCA <- as.princomp(phyl.pca(xipho_tree,xipho_phen_m_mat))
# as.princomp does not create a n.obs variable which is required for ggbiplot
xipho_phen_m_pPCA$n.obs <- length(rownames(xipho_phen_m_mat))

# write.csv(xipho_phen_m_pPCA$loadings, "data/male_pPCA_loadings_20240515.csv", row.names = T)

# add phylomorphospace to PCA 
# adapted from Liam Revell http://blog.phytools.org/2022/12/generating-ggplot-phylomorphospace-plot.html

# extract PC1 and PC2 and store in dataframe for ggobject
phylomorpho_m_comp <- xipho_phen_m %>%
  ungroup() %>%
  mutate(PC1 = xipho_phen_m_pPCA[["scores"]][,1],
         PC2 = xipho_phen_m_pPCA[["scores"]][,2]) %>%
  dplyr::select(PC1,PC2)
rownames(phylomorpho_m_comp) <- xipho_phen_m$tree_lab
phylomorpho_m_comp <- as.matrix(phylomorpho_m_comp)

phylomorpho_m_comp_obj <- phylomorphospace(xipho_tree,
                                           phylomorpho_m_comp,
                                           node.size=c(0,1.2),
                                           ftype="off")

# get coordinates of tree edges
phylomorpho_m_coords <- data.frame(
  xstart=phylomorpho_m_comp_obj$xx[phylomorpho_m_comp_obj$edge[,1]],
  ystart=phylomorpho_m_comp_obj$yy[phylomorpho_m_comp_obj$edge[,1]],
  xstop=phylomorpho_m_comp_obj$xx[phylomorpho_m_comp_obj$edge[,2]],
  ystop=phylomorpho_m_comp_obj$yy[phylomorpho_m_comp_obj$edge[,2]],
  nodestart=phylomorpho_m_comp_obj$edge[,1],
  nodestop=phylomorpho_m_comp_obj$edge[,2])

p_m <- ggbiplot(xipho_phen_m_pPCA,
                    var.axes = F,
                    ellipse = F,
                    group=xipho_phen_m$tactic,
                    segment.colour = 'grey',
                    segment.alpha = 0.5,
                    scale = 0) +
  scale_color_manual(values = c("courtship" = "#71018c", "coercive" = "#01a85b", "ambiguous" = "darkgrey")) +
  geom_segment(data=phylomorpho_m_coords,aes(x=xstart,y=ystart,xend=xstop,yend=ystop), 
               size = 0.5,color="grey",alpha=0.5) +
  geom_point(aes(colour= xipho_phen_m$tactic), pch = 16, size=3) +
  theme_bw() + 
  ggtitle("Average male phenotypes") +
  labs(color = "male tactic") +
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_text_repel(label=xipho_phen_m$label, point.padding = 0.5, box.padding = 0.2)

###################
# female analysis #
###################
xipho_phen_f <- xipho_phen %>%
  ungroup() %>%
  filter(sex == 1) %>%
  dplyr::select(-c(sex,
                   swd_length,
                   vertical_bars,
                   caudal_fin_pigmentation,
                   dorsal_fin_pigmentation,
                   upper_sword_edge_width,
                   upper_sword_edge_length,
                   lower_sword_edge_width,
                   lower_sword_edge_length,
                   dorsal_fin_length,
                   dorsal_fin_height,
                   body_depth,
                   peduncle_depth,
                   caudal_fin_length,
                   caudal_fin_height,
                   false_gravid_spot,
                   peduncle_edge,
                   lower_edge_2_peduncle))
rownames(xipho_phen_f) <- xipho_phen_f$tree_lab

# collapse polymorphic males into one tip for female analysis
xipho_tree_f <- drop.tip(xipho_tree, c("Xmul_small","Xnig_small"))

xipho_phen_f <- xipho_phen_f %>%
  arrange(factor(tree_lab, levels = xipho_tree_f$tip.label)) %>%
  # remove sword traits and atro since there were no females with this pattern
  dplyr::select(-c(atro,
                   upper_sword_edge_width_residuals,
                   upper_sword_edge_length_residuals,
                   lower_sword_edge_length_residuals,
                   lower_sword_edge_width_residuals,
                   sword_presence,
                   sword_pigmentation,
                   swd_length_residuals,
                   lower_edge_2_peduncle_residuals)) %>%
  dplyr::mutate(label = case_when(tree_lab %in% c("Xnig_small","Xnig_large","Xmul_small","Xmul_large","Xcon","Xmon","Xpyg","Xbir","Xhel") ~ tree_lab,
                           !(tree_lab %in% c("Xnig_small","Xnig_large","Xmul_small","Xmul_large","Xcon","Xmon","Xpyg","Xbir","Xhel")) ~ "")) %>%
  dplyr::mutate(label = case_when(label == "Xmul_large" ~ "Xmul",
                                  label == "Xnig_large" ~ "Xnig",
                                  label != "Xmul_large" & label != "Xnig_large" ~ label))
rownames(xipho_phen_f) <- xipho_phen_f$tree_lab

# keep one instance of this data as a matrix with rownames for things in phytools
xipho_phen_f_mat <- as.matrix(xipho_phen_f[5:16])
rownames(xipho_phen_f_mat) <- rownames(xipho_phen_f)

# run phylogenetic PCA using female tree (no sneaker morphs)
xipho_phen_f_pPCA <- as.princomp(phyl.pca(xipho_tree_f,xipho_phen_f_mat))
# as.princomp does not create a n.obs variable which is required for ggbiplot
xipho_phen_f_pPCA$n.obs <- length(rownames(xipho_phen_f_mat))

# write.csv(xipho_phen_f_pPCA$loadings, "data/female_pPCA_loadings_20240515.csv", row.names = T)

# add phylomorphospace to PCA
# extract PC1 and PC2 and store in dataframe for ggobject
phylomorpho_f_comp <- xipho_phen_f %>%
  ungroup() %>%
  mutate(PC1 = xipho_phen_f_pPCA[["scores"]][,1],
         PC2 = xipho_phen_f_pPCA[["scores"]][,2]) %>%
  dplyr::select(PC1,PC2)
rownames(phylomorpho_f_comp) <- xipho_phen_f$tree_lab
phylomorpho_f_comp <- as.matrix(phylomorpho_f_comp)

phylomorpho_f_comp_obj <- phylomorphospace(xipho_tree_f,
                                           phylomorpho_f_comp,
                                           node.size=c(0,1.2),
                                           ftype="off")

# get coordinates of tree edges
phylomorpho_f_coords <- data.frame(
  xstart=phylomorpho_f_comp_obj$xx[phylomorpho_f_comp_obj$edge[,1]],
  ystart=phylomorpho_f_comp_obj$yy[phylomorpho_f_comp_obj$edge[,1]],
  xstop=phylomorpho_f_comp_obj$xx[phylomorpho_f_comp_obj$edge[,2]],
  ystop=phylomorpho_f_comp_obj$yy[phylomorpho_f_comp_obj$edge[,2]],
  nodestart=phylomorpho_f_comp_obj$edge[,1],
  nodestop=phylomorpho_f_comp_obj$edge[,2])

p_f <- ggbiplot(xipho_phen_f_pPCA,
                    var.axes = F,
                    ellipse = F,
                    group=xipho_phen_f$tactic,
                    segment.colour = 'grey',
                    segment.alpha = 0.5,
                    scale = 0) +
  geom_segment(data=phylomorpho_f_coords,aes(x=xstart,y=ystart,xend=xstop,yend=ystop), 
               size = 0.5,color="grey",alpha=0.5) +
  geom_point(aes(colour= xipho_phen_f$tactic), pch = 16, size=3) +
  scale_color_manual(values = c("courtship" = "#71018c", "coercive" = "#01a85b", "ambiguous" = "darkgrey", "polymorphic" = "#FCC000")) +
  theme_bw() + 
  ggtitle("Average female phenotypes") +
  labs(color = "male tactic") +
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_text_repel(label=xipho_phen_f$label, point.padding = 0.5, box.padding = 0.2)

##############################
# sexual dimorphism analysis #
##############################

# get species names to eventually cbind with dimorphism matrix
sex_di_species <- xipho_phen %>%
  ungroup() %>%
  filter(sex == 0) %>%
  dplyr::select(species,tactic,clade,tree_lab)

# get male trait matrix
sex_di_m <- xipho_phen %>%
  ungroup() %>%
  dplyr::filter(sex == 0) %>%
  dplyr::select(-c(species,sex,tactic,clade,tree_lab))

# get female trait matrix (and keep grouping vars for now)
sex_di_f <- xipho_phen %>%
  dplyr::filter(sex == 1)

# extract xmul and xnig females to duplicate in next step
sex_di_f_dup <- sex_di_f %>%
  filter(species == "Xmultilineatus" | species == "Xnigrensis")

# bind extra xmul and xnig rows, sort by species
# see that there are 2 duplicated rows for xmul and xnig
sex_di_f <- rbind(sex_di_f,sex_di_f_dup) %>%
  arrange(species) %>%
  ungroup() %>%
  dplyr::select(-c(species,sex,tactic,clade,tree_lab)) # get rid of grouping vars for m/f subtraction

# subtract male trait matrix from female trait matrix
sex_di_for_pca <- cbind(sex_di_species, (sex_di_m - sex_di_f))

# filter out non-size-corrected variables
sex_di_for_pca <- sex_di_for_pca %>%
  dplyr::select(-c(swd_length,
                   vertical_bars,
                   caudal_fin_pigmentation,
                   dorsal_fin_pigmentation,
                   upper_sword_edge_width,
                   upper_sword_edge_length,
                   lower_sword_edge_width,
                   lower_sword_edge_length,
                   dorsal_fin_length,
                   dorsal_fin_height,
                   body_depth,
                   peduncle_depth,
                   caudal_fin_length,
                   caudal_fin_height,
                   peduncle_edge,
                   lower_edge_2_peduncle))
rownames(sex_di_for_pca) <- sex_di_for_pca$tree_lab

# match rownames with tree tip order
sex_di_for_pca <- sex_di_for_pca %>%
  arrange(factor(tree_lab, levels = xipho_tree$tip.label)) %>%
  # add label to only highlight northern swordtails
  mutate(label = case_when(tree_lab %in% c("Xnig_small","Xnig_large","Xmul_small","Xmul_large","Xcon","Xmon","Xpyg","Xhel","Xbir") ~ tree_lab,
                           !(tree_lab %in% c("Xnig_small","Xnig_large","Xmul_small","Xmul_large","Xcon","Xmon","Xpyg","Xhel","Xbir")) ~ "")) %>%
  mutate(label=str_replace(label,"_.*",""))
rownames(sex_di_for_pca) <- sex_di_for_pca$tree_lab

# make sex di matrix class obj
sex_di_for_pca_mat <- as.matrix(sex_di_for_pca[5:26])
rownames(sex_di_for_pca_mat) <- rownames(sex_di_for_pca)

# run phylogenetic PCA on male-female diff
sex_di_pPCA <- as.princomp(phyl.pca(xipho_tree,sex_di_for_pca_mat))
# as.princomp does not create a n.obs variable which is required for ggbiplot
sex_di_pPCA$n.obs <- length(rownames(sex_di_for_pca))

# write.csv(sex_di_pPCA$loadings, "data/sex_di_pPCA_loadings_20240515.csv", row.names = T)

# add phylomorphospace to PCA
# extract PC1 and PC2 and store in dataframe for ggobject
phylomorpho_sex_di_comp <- sex_di_for_pca %>%
  ungroup() %>%
  mutate(PC1 = sex_di_pPCA[["scores"]][,1],
         PC2 = sex_di_pPCA[["scores"]][,2]) %>%
  dplyr::select(PC1,PC2)
rownames(phylomorpho_sex_di_comp) <- sex_di_for_pca$tree_lab
phylomorpho_sex_di_comp <- as.matrix(phylomorpho_sex_di_comp)

phylomorpho_sex_di_comp_obj <- phylomorphospace(xipho_tree,
                                                phylomorpho_sex_di_comp,
                                                node.size=c(0,1.2),
                                                ftype="off")

# get coordinates of tree edges
phylomorpho_sex_di_coords <- data.frame(
  xstart=phylomorpho_sex_di_comp_obj$xx[phylomorpho_sex_di_comp_obj$edge[,1]],
  ystart=phylomorpho_sex_di_comp_obj$yy[phylomorpho_sex_di_comp_obj$edge[,1]],
  xstop=phylomorpho_sex_di_comp_obj$xx[phylomorpho_sex_di_comp_obj$edge[,2]],
  ystop=phylomorpho_sex_di_comp_obj$yy[phylomorpho_sex_di_comp_obj$edge[,2]],
  nodestart=phylomorpho_sex_di_comp_obj$edge[,1],
  nodestop=phylomorpho_sex_di_comp_obj$edge[,2])

# make extra point to plot on top of X. continens point
# messing with the factor order within the princomp object was giving me trouble so this is my workaround
sex_di_Xcon <- data.frame(Comp.1=c(sex_di_pPCA$scores[11,1]),
                          Comp.2=c(sex_di_pPCA$scores[11,2]))

# plot sexual dimorphism PCA coded by tactic
p_sex_di <- ggbiplot(sex_di_pPCA,
                         var.axes = F,
                         ellipse = F,
                         group=sex_di_for_pca$tactic,
                         segment.colour = 'grey',
                         segment.alpha = 0.5,
                         scale=0) +
  geom_segment(data=phylomorpho_sex_di_coords,aes(x=xstart,y=ystart,xend=xstop,yend=ystop),
  size = 0.5,color="grey",alpha=0.5) +
  scale_color_manual(values = c("courtship" = "#71018c", "coercive" = "#01a85b", "ambiguous" = "darkgrey")) +
  geom_point(aes(colour= sex_di_for_pca$tactic), pch = 16, size=3) +
  geom_point(aes(x = sex_di_Xcon$Comp.1, y= sex_di_Xcon$Comp.2), pch = 16, size=3, color="#01a85b") +
  # need to do this piecewise bc "! Aesthetics must be either length 1 or the same as the data (26)"
  ggtitle("Sexual dimorphism") +
  labs(color = "male tactic") +
  theme_bw() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_text_repel(label=sex_di_for_pca$label, point.padding = 0.5, box.padding = 0.2)

##########################################################
# PGLS between standard length and sexual dimorphism PC1 #
##########################################################

# add PC1 to male data
xipho_phen_m_comp <- xipho_phen_m %>%
  mutate(PC1 = sex_di_pPCA[["scores"]][,1],
         PC2 = sex_di_pPCA[["scores"]][,2])
rownames(xipho_phen_m_comp) <- xipho_phen_m_comp$tree_lab
#write.csv(xipho_phen_m_comp, file="data/xipho_phen_m_comp_PC1-2.csv")

xipho_phen_m_comp <- comparative.data(xipho_tree,
                                      as.data.frame(xipho_phen_m_comp),
                                      names.col = 'tree_lab',
                                      vcv.dim = 2,
                                      warn.dropped = T)

m_pgls_std_length_v_sex_di_PC1_caper <- pgls(PC1~std_length, data=xipho_phen_m_comp)
summary(m_pgls_std_length_v_sex_di_PC1_caper)
# significant ***

#####################################
# within-species variation male PCA #
#####################################

#perform similar filtering as xipho_phen but don't average
indiv_m_phen <- indiv_m_phen %>%
  dplyr::filter(!species %in% c("Xmonticolus", "Pjonesii")) %>%
  dplyr::filter(sex == "M") %>%
  mutate(sword_index = swd_length/std_length,
         # compute sword presence variable
         sword_presence = case_when(sword_index > 0.01943906 ~ 1,
                                    sword_index <= 0.01943906 ~ 0)) %>%
  dplyr::select(-c(sex,
                   fish_ID,
                   swd_length,
                   vertical_bars,
                   caudal_fin_pigmentation,
                   dorsal_fin_pigmentation,
                   upper_sword_edge_width,
                   upper_sword_edge_length,
                   lower_sword_edge_width,
                   lower_sword_edge_length,
                   dorsal_fin_length,
                   dorsal_fin_height,
                   body_depth,
                   peduncle_depth,
                   caudal_fin_length,
                   caudal_fin_height,
                   peduncle_edge,
                   lower_edge_2_peduncle))

indiv_m_phen_NS <- indiv_m_phen %>%
  filter(clade=="NS" & tree_lab %in% c("Xcon","Xpyg","Xmon","Xmul_large","Xmul_small","Xnig_large","Xnig_small")) %>%
  mutate(tree_lab = case_when(tree_lab %in% c("Xmul_large","Xmul_small") ~ "Xmul",
                             tree_lab %in% c("Xnig_large","Xnig_small") ~ "Xnig",
                             !tree_lab %in% c("Xnig_large","Xnig_small","Xmul_large","Xmul_small") ~ tree_lab))
indiv_m_phen_NS_mat <- as.matrix(indiv_m_phen_NS %>%
                                   dplyr::select(-c(species,tactic,clade,sword_index,tree_lab,atro,false_gravid_spot,sc)))

indiv_m_NS_pca <- princomp(indiv_m_phen_NS_mat)
summary(indiv_m_NS_pca)

# write.csv(indiv_m_NS_pca$loadings, "data/indiv_m_NS_pca_loadings_20240515.csv", row.names = T)

# doing this with ggplot -- using shapes as an additional aesthetic in ggbiplot was giving me trouble
indiv_m_NS_pca_dat <- as.data.frame(indiv_m_NS_pca$scores) %>%
  dplyr::select(Comp.1,Comp.2)
indiv_m_NS_pca_dat <- cbind(indiv_m_phen_NS, indiv_m_NS_pca_dat)
indiv_m_NS_pca_dat <- indiv_m_NS_pca_dat %>%
  mutate(tactic = case_when(tactic == "sneaker" ~ "coercive",
                            tactic == "courter" ~ "courtship",
                            tactic %in% c("ambiguous","polymorphic") ~ tactic))

p_ins <- ggplot(indiv_m_NS_pca_dat, aes(x=Comp.1,y=Comp.2,color=as.factor(tactic))) +
  geom_point(aes(shape=as.factor(species)), size=3) +
  scale_color_manual(values = c("courtship" = "#71018c", "coercive" = "#01a85b", "ambiguous" = "darkgrey")) +
  theme_bw() +
  labs(color="male tactic",
       shape="species",
       x="PC1 (86.4% explained var.)",
       y="PC2 (10.1% explained var.)") +
  ggtitle("Within-species male phenotypic variation") +
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

################################
# SL distribution joy div plot #
################################
sl_dists <- read.csv("data/sl_dists.csv")

p_sl <- ggplot(sl_dists, aes(x=std_length, y=species)) +
  geom_density_ridges(alpha=0.8) +
  labs(x="standard length (mm)") +
  ggtitle("Population level male size distribution") +
  theme_test() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()) +
  annotate("text", x = 13,y = 5+0.5,label = "Xvar ", fontface="italic") +
  annotate("text", x = 13,y = 4+0.5,label = "Xnig ",fontface="italic") +
  annotate("text", x = 13,y = 3+0.5,label = "Xnez ",fontface="italic") +
  annotate("text", x = 13,y = 2+0.5,label = "Xmul *",fontface="italic") +
  annotate("text", x = 13,y = 1+0.5,label = "Xcor ",fontface="italic")

dip.test((sl_dists %>% filter(species == "X. variatus"))$std_length)
dip.test((sl_dists %>% filter(species == "X. nigrensis"))$std_length)
dip.test((sl_dists %>% filter(species == "X. nezahualcoyotl"))$std_length)
dip.test((sl_dists %>% filter(species == "X. multilineatus"))$std_length)
dip.test((sl_dists %>% filter(species == "X. cortezi"))$std_length)

#########################################################
# make body size loss cladogram for northern swordtails #
#########################################################

# make separate tree to modify for visualization purposes
xipho_tree_cp <- as_tibble(xipho_tree) %>%
  mutate(label = case_when(label == "Xmal" ~ "X. malinche",
                           label=="Xbir" ~ "X. birchmanni",
                           label=="Xcor" ~ "X. cortezi",
                           label=="Xmon" ~ "X. montezumae",
                           label=="Xcon" ~ "X. continens",
                           label=="Xnez" ~ "X. nezahualcoyotl",
                           label=="Xmul_small" ~ "X. multilineatus",
                           label=="Xmul_large" ~ "X. multilineatus",
                           label=="Xnig_small" ~ "X. nigrensis",
                           label=="Xnig_large" ~ "X. nigrensis",
                           label=="Xpyg" ~ "X. pygmaeus",
                           label=="Xmey" ~ "X. meyeri",
                           label=="Xgor" ~ "X. gordoni",
                           label=="Xcou" ~ "X. couchianus",
                           label=="Xvar" ~ "X. variatus",
                           label=="Xeve" ~ "X. evelynae",
                           label=="Xmil" ~ "X. milleri",
                           label=="Xxip" ~ "X. xiphidium",
                           label=="Xand" ~ "X. andersi",
                           label=="Xmac" ~ "X. maculatus",
                           label=="Xsig" ~ "X. signum",
                           label=="Xalv" ~ "X. alvarezi",
                           label=="Xmay" ~ "X. mayae",
                           label=="Xkal" ~ "X. kallmani",
                           label=="Xhel" ~ "X. hellerii",
                           label=="Xcle" ~ "X. clemenciae",
                           label=="Gamb" ~ "G. affinis"))
xipho_tree_cp <- as.phylo(xipho_tree_cp)

p_tree <- ggtree(xipho_tree_cp, branch.length = "none", size=1) +
  geom_tiplab(hjust=-0.1, fontface="italic") + 
  #geom_text(aes(label=node)) +
  geom_point2(aes(subset=(label %in% c("X. continens","X. pygmaeus"))), color ="#01a85b", shape = 15, size=4) +
  geom_point2(aes(subset=(label == "X. multilineatus" & node == "9")), color ="#01a85b", shape = 15, size=4) +
  geom_point2(aes(subset=(label == "X. multilineatus" & node == "10")), color ="#71018C", shape = 15, size=4) +
  geom_point2(aes(subset=(label == "X. nigrensis" & node == "7")), color ="#01a85b", shape = 15, size=4) +
  geom_point2(aes(subset=(label == "X. nigrensis" & node == "8")), color ="#71018C", shape = 15, size=4) +
  geom_point2(aes(subset=(label %in% c("X. montezumae", "X. malinche", "X. cortezi","X. nezahualcoyotl", "X. birchmanni"))), color ="#71018C", shape = 15, size=4) +
  xlim(0,35) +
  theme(aspect.ratio = 1)
p_tree <- p_tree %>%
  scaleClade(node=44, 0.2) %>%
  scaleClade(node=29, 0.2) %>%
  ggtree::collapse(node=44, "max", color="black", fill="white", size=1) %>%
  ggtree::collapse(node=29, "max", color="black", fill="white", size=1)
p_tree <- p_tree +
  geom_cladelabel(node=44, "platyfish", hjust=-1.8, vjust=-1.5) +
  geom_cladelabel(node=29, "southern swordtails", hjust=-.4, vjust=-.25)

##############################
# cowplot everything togther #
##############################

#Fig 2
p_tree_sex_di <- plot_grid(p_tree, p_sex_di, labels=c("A","B"), nrow=1)
p_ins_sl <- plot_grid(p_ins, p_sl, labels=c("C","D"), nrow=1)

plot_grid(
  p_tree_sex_di, p_ins_sl,
  nrow=2
)

# Fig S6
plot_grid(
  p_m, p_f,
  nrow=1, labels="AUTO"
)