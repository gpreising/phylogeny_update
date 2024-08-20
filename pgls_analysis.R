# This script creates Fig 3 and shows the general process we used for simulations in Fig S10-12
# Please cite "Recurrent evolution of small body size and loss of the sword ornament in Northern swordtail fish" Preising et al. 2024  Evolution if you use this script

library(tidyverse)
library(ape)
library(phytools)
library(caper)
library(cowplot)

# read in data
xipho_phen <- read.csv("data/phy_test_info_20240430_recode.csv") %>%
  filter(sex=="M")
# make vector of allometric traits to iterate over
res_traits <- c("swd_length","upper_sword_edge_width","upper_sword_edge_length","lower_sword_edge_width","lower_sword_edge_length",
                "dorsal_fin_length","dorsal_fin_height","body_depth","peduncle_depth","caudal_fin_height","caudal_fin_length",
                "peduncle_edge","lower_edge_2_peduncle","vertical_bars")

# compute residuals for allometric traits and store as vars in initial df
for (trait in res_traits){
  
  formula <- as.formula(paste(trait,"~ std_length"))
  lm_model <- lm(formula,data=xipho_phen)
  residuals <- lm_model$residuals
  varname <- paste(trait, "_residuals", sep = "")
  
  xipho_phen <- xipho_phen %>%
    dplyr::mutate(!!varname := residuals)
  
}

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
  dplyr::select(c(species,
                  sex,
                  std_length,
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
                  vertical_bars_residuals))

#####################
# read in tree data #
#####################

# read in tree
xipho_tree <- read.nexus("data/RAxML_bipartitions.NS_complete_20231215_treefile_addmorphs.txt.tre")
# remove kallmani
xipho_tree <- drop.tip(xipho_tree, c("Xkal"))

###############################
# prep data for PGLS analysis #
###############################

# create a subset for males and make the names in the dataframe match the names in the tree
xipho_phen_m <- as.data.frame(xipho_phen %>%
  filter(sex == 0) %>%
  dplyr::mutate(
    species = case_when(species == "Gaffinis" ~ "Gamb",
                        species == "Xalvarezi" ~ "Xalv",
                        species == "Xandersi" ~ "Xand",
                        species == "Xbirchmanni" ~ "Xbir",
                        species == "Xclemenciae" ~ "Xcle",
                        species == "Xcontinens" ~ "Xcon",
                        species == "Xcortezi" ~ "Xcor",
                        species == "Xcouchianus" ~ "Xcou",
                        species == "Xevelynae" ~ "Xeve",
                        species == "Xgordoni" ~ "Xgor",
                        species == "Xhellerii" ~ "Xhel",
                        species == "Xmaculatus" ~ "Xmac",
                        species == "Xmalinche" ~ "Xmal",
                        species == "Xmayae" ~ "Xmay",
                        species == "Xmeyeri" ~ "Xmey",
                        species == "Xmilleri" ~ "Xmil",
                        species == "Xmontezumae" ~ "Xmon",
                        species == "Xmultilineatus" & tactic == "courter" ~ "Xmul_large",
                        species == "Xmultilineatus" & tactic == "sneaker" ~ "Xmul_small",
                        species == "Xnezahualcoyotl" ~ "Xnez",
                        species == "Xnigrensis" & tactic == "courter" ~ "Xnig_large",
                        species == "Xnigrensis" & tactic == "sneaker" ~ "Xnig_small",
                        species == "Xpygmaeus" ~ "Xpyg",
                        species == "Xsignum" ~ "Xsig",
                        species == "Xvariatus" ~ "Xvar",
                        species == "Xxiphidium" ~ "Xxip")) %>%
  ungroup() %>%
  dplyr::select(-c(sc,sex,atro,tactic,clade)) %>%
  arrange(factor(species, levels = xipho_tree$tip.label)))
rownames(xipho_phen_m) <- xipho_phen_m$species

################
# PGLS - caper # 
################

# create comparative.data object for caper
xipho_comp_data <- comparative.data(xipho_tree,
                                    xipho_phen_m,
                                    names.col = 'species',
                                    vcv.dim = 2,
                                    warn.dropped = T)

# sword length residuals
m_pgls_std_length_v_swd_length_caper <- pgls(swd_length_residuals~std_length, data=xipho_comp_data)
summary(m_pgls_std_length_v_swd_length_caper)

# upper sword edge width residuals
m_pgls_upperedgewidth_v_sl_length_caper <- pgls(upper_sword_edge_width_residuals~std_length, data=xipho_comp_data)
summary(m_pgls_upperedgewidth_v_sl_length_caper)

# upper sword edge length residuals
m_pgls_upperedgelength_sl_length_caper <- pgls(upper_sword_edge_length_residuals~std_length, data=xipho_comp_data)
summary(m_pgls_upperedgelength_sl_length_caper)

# dorsal fin length residuals
m_pgls_dorsallength_v_sl_length_caper <- pgls(dorsal_fin_length_residuals~std_length, data=xipho_comp_data)
summary(m_pgls_dorsallength_v_sl_length_caper)

# sword presence - marginally insignificant .
m_pgls_sp_v_sl_length_caper <- pgls(sword_presence ~ std_length, data=xipho_comp_data)
summary(m_pgls_sp_v_sl_length_caper)

# caudal fin pigmentation
m_pgls_cfp_sl_length_caper <- pgls(caudal_fin_pigmentation~std_length, data=xipho_comp_data)
summary(m_pgls_cfp_sl_length_caper)

# lower sword edge width residuals
m_pgls_lsew_sl_length_caper <- pgls(lower_sword_edge_width_residuals~std_length, data=xipho_comp_data)
summary(m_pgls_lsew_sl_length_caper)

# lower sword edge length residuals
m_pgls_lsel_sl_length_caper <- pgls(lower_sword_edge_length_residuals~std_length, data=xipho_comp_data)
summary(m_pgls_lsel_sl_length_caper)

# caudal blotch - margnially insignificant .
m_pgls_cb_sl_length_caper <- pgls(cb~std_length, data=xipho_comp_data)
summary(m_pgls_cb_sl_length_caper)

# false gravid spot
m_pgls_fgs_sl_length_caper <- pgls(false_gravid_spot~std_length, data=xipho_comp_data)
summary(m_pgls_fgs_sl_length_caper)

# carbomaculatus - marginally insignificant .
m_pgls_carbo_sl_length_caper <- pgls(carbo~std_length, data=xipho_comp_data)
summary(m_pgls_carbo_sl_length_caper)

# dorsal fin pigmentation
m_pgls_dpig_v_sl_length_caper <- pgls(dorsal_fin_pigmentation ~std_length, data=xipho_comp_data)
summary(m_pgls_dpig_v_sl_length_caper)

# caudal fin length residuals
m_pgls_cfl_sl_length_caper <- pgls(caudal_fin_length_residuals~std_length, data = xipho_comp_data)
summary(m_pgls_cfl_sl_length_caper)

# peduncle edge residuals - significant *
m_pgls_pe_sl_length_caper <- pgls(peduncle_edge_residuals~std_length, data=xipho_comp_data)
summary(m_pgls_pe_sl_length_caper)

# caudal fin height residuals - significant *
m_pgls_cfh_sl_length_caper <- pgls(caudal_fin_height_residuals~std_length, data=xipho_comp_data)
summary(m_pgls_cfh_sl_length_caper)

# sword pigmentation - significant *
m_pgls_spig_v_sl_length_caper <- pgls(sword_pigmentation~std_length, data=xipho_comp_data)
summary(m_pgls_spig_v_sl_length_caper)

# vertical bars - non body size corrected - significant ***
m_pgls_vb_v_sl_length_caper <- pgls(vertical_bars~std_length, data=xipho_comp_data)
summary(m_pgls_vb_v_sl_length_caper)

# dorsal fin height residuals - significant ***
m_pgls_dorsalheight_v_sl_length_caper <- pgls(dorsal_fin_height_residuals~std_length, data=xipho_comp_data)
summary(m_pgls_dorsalheight_v_sl_length_caper)
p_dfh <- ggplot(xipho_phen_m, aes(x = std_length, y=dorsal_fin_height_residuals)) +
  geom_abline(intercept = m_pgls_dorsalheight_v_sl_length_caper[["model"]][["coef"]][[1]],
              slope = m_pgls_dorsalheight_v_sl_length_caper[["model"]][["coef"]][[2]],
              col="gray60",
              lty=2,
              lwd=1) +
  geom_point(size=4,alpha=0.6,col="blue") + 
  theme_bw() +
  labs(x = "standard length (mm)",
       y = "dorsal fin height residuals") +
  theme(aspect.ratio = 1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  annotate("text", x = min(xipho_phen_m$std_length)+6,
           y = min(xipho_phen_m$dorsal_fin_height_residuals),
           label = paste("P = ",round(summary(m_pgls_dorsalheight_v_sl_length_caper)[["coefficients"]][[8]], 4)))

# body depth residuals - significant ***
m_pgls_bodydepth_v_sl_length_caper <- pgls(body_depth_residuals~std_length, data=xipho_comp_data)
summary(m_pgls_bodydepth_v_sl_length_caper)
p_bd <- ggplot(xipho_phen_m, aes(x = std_length, y=body_depth_residuals)) +
  geom_abline(intercept = m_pgls_bodydepth_v_sl_length_caper[["model"]][["coef"]][[1]],
              slope = m_pgls_bodydepth_v_sl_length_caper[["model"]][["coef"]][[2]],
              col="gray60",
              lty=2,
              lwd=1) +
  geom_point(size=4,alpha=0.6,col="blue") + 
  theme_bw() +
  labs(x = "standard length (mm)",
       y = "body depth residuals") +
  theme(aspect.ratio = 1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  annotate("text", x = min(xipho_phen_m$std_length)+6,
           y = min(xipho_phen_m$body_depth_residuals),
           label = paste("P = ",round(summary(m_pgls_bodydepth_v_sl_length_caper)[["coefficients"]][[8]], 4)))

# peduncle depth residuals - significant ***
m_pgls_pd_sl_length_caper <- pgls(peduncle_depth_residuals~std_length, data = xipho_comp_data)
summary(m_pgls_pd_sl_length_caper)
p_pd <- ggplot(xipho_phen_m, aes(x = std_length, y=peduncle_depth_residuals)) +
  geom_abline(intercept = m_pgls_pd_sl_length_caper[["model"]][["coef"]][[1]],
              slope = m_pgls_pd_sl_length_caper[["model"]][["coef"]][[2]],
              col="gray60",
              lty=2,
              lwd=1) +
  geom_point(size=4,alpha=0.6,col="blue") + 
  theme_bw() +
  labs(x = "standard length (mm)",
       y = "peduncle depth residuals") +
  theme(aspect.ratio = 1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  annotate("text", x = min(xipho_phen_m$std_length)+6,
           y = min(xipho_phen_m$peduncle_depth_residuals),
           label = paste("P = ",round(summary(m_pgls_pd_sl_length_caper)[["coefficients"]][[8]], 5)))

# vertical bars residuals - significant **
m_pgls_vb_res_sl_length_caper <- pgls(vertical_bars_residuals~std_length, data=xipho_comp_data)
summary(m_pgls_vb_res_sl_length_caper)
#plot(m_pgls_vb_res_sl_length_caper)
p_vb_res <- ggplot(xipho_phen_m, aes(x = std_length, y=vertical_bars_residuals)) +
  geom_abline(intercept = m_pgls_vb_res_sl_length_caper[["model"]][["coef"]][[1]],
              slope = m_pgls_vb_res_sl_length_caper[["model"]][["coef"]][[2]],
              col="gray60",
              lty=2,
              lwd=1) +
  geom_point(size=4,alpha=0.6,col="blue") + 
  theme_bw() +
  labs(x = "standard length (mm)",
       y = "vertical bars residuals") +
  theme(aspect.ratio = 1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  annotate("text", x = min(xipho_phen_m$std_length)+6,
           y = min(xipho_phen_m$vertical_bars_residuals),
           label = paste("P = ",round(summary(m_pgls_vb_res_sl_length_caper)[["coefficients"]][[8]], 4)))

# read in pPCA data
xipho_phen_m_pPC1_2 <- read.csv("data/xipho_phen_m_comp_PC1-2.csv")
xipho_phen_m_pPC1_2_comp <- comparative.data(xipho_tree,
                                             xipho_phen_m_pPC1_2,
                                             names.col = 'species',
                                             vcv.dim = 2,
                                             warn.dropped = T)
# pPC1 - significant ***
m_pgls_pPC1_sl_length_caper <- pgls(PC1~std_length, data=xipho_phen_m_pPC1_2_comp)
summary(m_pgls_pPC1_sl_length_caper)

###########
# cowplot #
###########

plot_grid(p_pd, p_bd, p_dfh, p_vb_res, nrow=2, ncol=2, labels='AUTO')


#########################
###########power simulations for continous traits
####BM mod with different variances
####drawn from here:http://phytools.org/eqg/Exercise_4.1/
#########################

#set gamma
set_gamma=0.2

#set l1 and l2

mod2_BM <- function(x, l1,l2, gamma) {
  l1=10
  l2=0.01
  out1 <- rnorm(1, x[1] + gamma * x[2], sqrt(l1 * 1))
  out2 <- rnorm(1, gamma * x[1] + x[2], sqrt(l2 * 1))
  c(out1, out2)
}


###run power simulations
r2<-{}
p<-{}

for(x in 1:1000){
  
  q<-rTraitMult(xipho_tree,mod2_BM,gamma=set_gamma,p=2,root.value=c(33,0)) #set root values trait by trait
  xipho_phen_m$sim1<-q[,1]
  xipho_phen_m$sim2<-q[,2]
  
  xipho_sim_data <- comparative.data(xipho_tree,
                                     xipho_phen_m,
                                     names.col = 'species',
                                     vcv.dim = 2,
                                     warn.dropped = T)
  
  sim_run <- pgls(sim1~sim2, data=xipho_sim_data)
  p<-c(p,summary(sim_run)$coefficients[8])
  r2<-c(r2,summary(sim_run)$r.squared)
  
}
length(subset(p,p<0.05))/1000
median(r2)

##########################################
###simulate continuous traits independently on a tree (false positive rate)
###########################################
x<-xipho_phen_m$std_length
#set alpha and sig2
alpha=0.025
sigtwo=15
y<-fastBM(xipho_tree, a=alpha, sig2=sigtwo,mu=0.1,bounds=c(-3,4)) 

#visualize one replicate versus real data (as an example, body depth)
par(mfrow=c(1,2))
hist(xipho_phen_m$body_depth_residuals,xlab="Observed residuals",main="",col="lightblue")
hist(y,xlab="Simulated residuals",main="",col="lightgray")

###run null simulations

r2<-{}
null_pvals<-{}

for(k in 1:1000){
  y<-fastBM(xipho_tree, a=alpha, sig2=sigtwo,mu=0.1,bounds=c(-3,4)) #body depth
  xipho_phen_m$sim1<-y
  xipho_sim_data <- comparative.data(xipho_tree,
                                     xipho_phen_m,
                                     names.col = 'species',
                                     vcv.dim = 2,
                                     warn.dropped = T)
  sim_run <- pgls(std_length~sim1, data=xipho_sim_data)
  null_pvals<-c(null_pvals,summary(sim_run)$coefficients[8])
  r2<-c(r2,summary(sim_run)$r.squared)
}

#check whether p-values are well-calibrated
length(subset(null_pvals,null_pvals<0.05))/1000

hist(null_pvals,xlab="Simulation p-value",ylab="Frequency",col="gray",main="Body depth simulation")

