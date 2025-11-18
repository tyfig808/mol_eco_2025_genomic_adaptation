## plot pca from snps in clear and cmh tests
setwd("~/mol_eco_2024/2.5.4_pca")

library(tidyverse)
library(vegan)
library(ggplot2)
if (!require(ggiraph)) install.packages('ggiraph')
library(ggpubr)
library(data.table)
if (!require(svglite)) install.packages('svglite')

### read in eigenvect from plink, load full first
full_pca <- fread(paste0("pca_final_total_snp_no_indel.eigenvec")) 
sig_pca <- fread(paste0("pca_clear_cmh_sig_no_indel.eigenvec"))


pca_list <- list(full_pca, sig_pca)
pca_names <- c("full_pca", "sig_pca")

pca_list <- lapply(pca_list, function(dt) {
  setnames(dt, "#IID", "id", skip_absent=TRUE)

  ### replace 27A with resequence 
  dt <- dt %>%
    #subset(sigzone_pca$V1 != "RA_TL_H_B_27A") %>%
    mutate(id = dplyr::recode(id, "RA_TL_H_B_27A_c" = "RA_TL_H_B_27A"))

    # get rep, then soil, remove the second soil if it is not G1
  dt$rep<- str_split_i(dt$id, "_", 1) 
  dt$soil<- str_split_i(dt$id, "_", 2)

  g_1 <- which(dt$soil!="G1")
  str_sub(dt$soil[g_1],2,2) <- ""  # remove the second L from all groups but G1

  dt$herb<- ifelse(dt$soil != "G1", str_split_i(dt$id, "_", 3), NA)
  dt$bee<- ifelse(dt$soil != "G1", str_split_i(dt$id, "_", 4), NA)

  # add treatment and group for plotting later
  dt[, treatment := paste(soil, herb, bee, sep = '_')]
  dt$treatment<- ifelse(dt$treatment == "G1_NA_NA", "G1", dt$treatment)

  dt[, group := paste(rep, treatment, sep = '_')]


  return(dt)
})

# convert back to dt
for(i in 1:(length(pca_names))) {
  assign(pca_names[i], pca_list[[i]])
}

### 
eg_full <- fread(paste0("pca_final_total_snp_no_indel.eigenval"))
eg_sig <- fread(paste0("pca_clear_cmh_sig_no_indel.eigenval")) 


# get percentages to explain variance
es_full<-as.numeric(sum(eg_full))
e1_full<-as.numeric(round(eg_full[1]/es_full, 3)*100)
e2_full<-as.numeric(round(eg_full[2]/es_full, 3)*100)

es_sig<-as.numeric(sum(eg_sig))
e1_sig<-as.numeric(round(eg_sig[1]/es_sig, 3)*100)
e2_sig<-as.numeric(round(eg_sig[2]/es_sig, 3)*100)


# add colors for plotting
# soil_colors <- c(Limestone = "#0071c1ff", Tuff = "#538234ff")

group_colors <- c(RA_G1 = "black", RA_L_H_B = "#0071c1ff", RA_L_H_H = "#0071c1ff", RA_L_NH_B = "#0071c1ff", RA_L_NH_H = "#0071c1ff",
                  RA_T_H_B = "#538234ff", RA_T_H_H = "#538234ff",RA_T_NH_B = "#538234ff", RA_T_NH_H = "#538234ff",
                  RB_G1 = "black", RB_L_H_B = "#0071c1ff", RB_L_H_H = "#0071c1ff",RB_L_NH_B = "#0071c1ff", RB_L_NH_H = "#0071c1ff",
                  RB_T_H_B = "#538234ff", RB_T_H_H = "#538234ff",RB_T_NH_B = "#538234ff",RB_T_NH_H = "#538234ff")
                     

group_shapes <- c(RA_G1 = 22, RA_L_H_B = 24, RA_L_H_H = 24,RA_L_NH_B = 21, RA_L_NH_H = 21, 
                  RA_T_H_B = 24, RA_T_H_H = 24,RA_T_NH_B = 21,RA_T_NH_H = 21,
                  RB_G1 = 22, RB_L_H_B = 24, RB_L_H_H = 24, RB_L_NH_B = 21, RB_L_NH_H = 21, 
                  RB_T_H_B = 24, RB_T_H_H = 24, RB_T_NH_B = 21, RB_T_NH_H = 21)

# get centroids for plotting
cen_full <- aggregate(cbind(full_pca$PC1,full_pca$PC2)~group,full_pca,mean)
cen_sig <- aggregate(cbind(sig_pca$PC1,sig_pca$PC2)~group,sig_pca,mean)


# get min and max for axes
p1_max <- max(full_pca$PC1, sig_pca$PC1)
p1_min <- min(full_pca$PC1, sig_pca$PC1)
p2_max <- max(full_pca$PC2, sig_pca$PC2)
p2_min <- min(full_pca$PC2, sig_pca$PC2)

# create plots sep then ggarrange together
p1 <- ggarrange(
  ggplot(full_pca, aes(x = PC1, y = PC2, shape = group, col = group, label=group)) +
    geom_point(size = 2) + 
    labs(x= paste0("PC1 (", e1_full, "% Variation Explained)"), y= paste0("PC2 (", e2_full, "% Variation Explained)")) +
    theme_bw(base_size = 19) +
    #coord_cartesian(xlim = c(p1_min, p1_max), ylim=c(p2_min, p2_max)) +
    theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1)) +
    scale_shape_manual(values=group_shapes) +
    scale_color_manual(values=group_colors) +
    geom_point(data = cen_full, aes(x = V1, y = V2), fill = group_colors, size=7, colour = "black"), #+
    #geom_text(hjust=0, vjust=0),
  ncol= 1, nrow = 1, common.legend = TRUE, legend="none")

p2 <- ggarrange(
  ggplot(sig_pca, aes(x = PC1, y = PC2, shape = group, col = group)) +
    geom_point(size = 2) + 
    labs(x= paste0("PC1 (", e1_sig, "% Variation Explained)"), y= paste0("PC2 (", e2_sig, "% Variation Explained)")) +
    theme_bw(base_size = 19) +
    #coord_cartesian(xlim = c(p1_min, p1_max), ylim=c(p2_min, p2_max)) +
    theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1)) +
    scale_shape_manual(values=group_shapes) +
    scale_color_manual(values=group_colors) +
    geom_point(data = cen_sig, aes(x = V1, y = V2), fill = group_colors, size=8, colour = "black"),
  ncol= 1, nrow = 1, common.legend = TRUE, legend="none")


b_h_pca <- ggarrange(p1, p2, ncol = 1, nrow = 2, labels = c("A","B"))

#b_h_pca <- annotate_figure(b_h_pca, 
  #top = text_grob("PCA of Treatment Combinations", color = "black", face = "bold", size = 19),
  #left = text_grob(paste0("PC2 (", e2, "%)"), color = "black", face = "bold", size = 19, rot = 90),
  #bottom = text_grob(paste0("PC1 (", e1, "%)"), color = "black", face = "bold", size = 19))

#ggsave(filename = "all_pca_sep.png", plot = b_h_pca, width = 40, height = 23.30, dpi = 1200, units = "cm")


# for full svg
resolution <- 144
svglite(paste0("pca_clear_cmh_upd_col.svg", sep=""), width = 1440/resolution, height = 1440/resolution)

p2

dev.off()

# for help labelling centroid in inkscape
ggarrange(
  ggplot(data = cen_sig, aes(x = V1, y = V2, shape = group, col = group, label = group)) +
    #geom_point(size = 2) + 
    labs(x= paste0("PC1 (", e1_sig, "% Variation Explained)"), y= paste0("PC2 (", e2_sig, "% Variation Explained)")) +
    theme_bw(base_size = 19) +
    coord_cartesian(xlim = c(p1_min, p1_max), ylim=c(p2_min, p2_max)) +
    theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1)) +
    scale_shape_manual(values=group_shapes) +
    scale_color_manual(values=group_colors) +
    geom_point(data = cen_sig, aes(x = V1, y = V2), fill = group_colors, size=7, colour = "black") +
    geom_text(hjust=0, vjust=0),
  ncol= 1, nrow = 1, common.legend = TRUE, legend="none")