if (!require(tidyverse)) install.packages('tidyverse')
if (!require(data.table)) install.packages('data.table')
if (!require(ggpubr)) install.packages('ggpubr')
if (!require(viridis)) install.packages('viridis')
if (!require(scales)) install.packages('scales')
if (!require(ggnewscale)) install.packages('ggnewscale')
if (!require(svglite)) install.packages('svglite')
if (!require(igraph)) install.packages('igraph')
if (!require(FNN)) install.packages('FNN')

# read in fst long
setwd("/shares/schiestl.systbot.uzh/variants/data/tfig/ch_1_sub/2.5.1_fst_pair")
full_fst <- fread(sep = "\t", paste0("total_rep_sep_fst_full_data.tsv"))

# get average across rep
fst <- full_fst[, mean(WEIGHTED_FST), by = c("pop_1", "pop_2")] 

# change name
setnames(fst, "V1", "mean_fst")


# convert to pairwise matrix
data2 <- data.frame(pop_1 = fst$pop_2, pop_2 = fst$pop_1, mean_fst = fst$mean_fst) # flip the pop 1 and two 
data3 <- as.data.frame(fst) 
data3<- data3[c("pop_1", "pop_2", "mean_fst")]
df2 <- rbind(data3, data2)
df3 <- as.data.frame.matrix(xtabs(mean_fst ~ ., df2))
rownames(df3) <- colnames(df3)

# run mds but for 2 dimension
mds <- cmdscale(df3, k = 2, eig = TRUE)

# get var explain by 1 and 2
sum_e <- sum(mds$eig)
e_1 <- mds$eig[1]/sum_e
e_2 <- mds$eig[2]/sum_e

# create data table for plotti9ng
mat <- data.table(cbind(x = mds$points[,1], y = mds$points[,2]))
mat$treat <- colnames(df3)

# make eigenvalues for plotting
ev_1 <- round(e_1*100, 0)
ev_2 <- round(e_2*100, 0)

# add color and shapes, same as PCA but group names are different
group_colors <- c(G1 = "black", LHB = "#0071c1ff", LHH = "#0071c1ff", LNHB = "#0071c1ff", LNHH = "#0071c1ff", 
		THB = "#538234ff", THH = "#538234ff",TNHB = "#538234ff", TNHH = "#538234ff")

group_shapes <- c(G1 = 22, LHB = 24, LHH = 24,LNHB = 21, LNHH = 21,
	THB = 24, THH = 24, TNHB = 21, TNHH = 21)

# plot and save
setwd("/shares/schiestl.systbot.uzh/variants/data/tfig/gwas")

resolution <- 144
svglite(paste0("fst_network_mds_upd_col.svg", sep=""), width = 1080/resolution, height = 1080/resolution)

ggarrange(
  ggplot(mat, aes(x = x, y = -y, shape = treat, col = treat)) +
    geom_point(fill = group_colors, size=10, colour = "black") + 
    labs(x= paste0("Genetic Distance (Fst) Axis 1 (", ev_1, "% Variation Explained)"), 
         y= paste0("Genetic Distance (Fst) Axis 2 (", ev_2, "% Variation Explained)")) +
    theme_bw(base_size = 19) +
    #coord_cartesian(xlim = c(p1_min, p1_max), ylim=c(p2_min, p2_max)) +
    theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1)) +
    scale_shape_manual(values=group_shapes) +
    scale_color_manual(values=group_colors),
  ncol= 1, nrow = 1, common.legend = TRUE, legend="none")


dev.off()

