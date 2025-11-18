# check the vif between phenotypic traits, how many would we have to drop
setwd("~/mol_eco_2024/2.4_pheno_evo")

if (!require(tidyverse)) install.packages('tidyverse')
if (!require(data.table)) install.packages('data.table')
if (!require(car)) install.packages('car')
if (!require(gammit)) install.packages('gammit')
if (!require(ggpubr)) install.packages('ggpubr')
if (!require(ggsignif)) install.packages('ggsignif')
if (!require(viridis)) install.packages('viridis')
if (!require(scales)) install.packages('scales')
if (!require(ggnewscale)) install.packages('ggnewscale')
if (!require(svglite)) install.packages('svglite')
if (!require(igraph)) install.packages('igraph')
if (!require(gridBase)) install.packages('gridBase')
if (!require(grid)) install.packages('grid')
if (!require(ggridges)) install.packages('ggridges')
# read in phen
phen <- fread(sep = "\t", paste0("emmax_pheno_data.tsv"))

# get correlation matrix, remove diag, and get which r above 0.7
x <- cor(na.omit(phen[,-1]) )
diag(x) <- NA
which(abs(x)>0.7, arr.ind = TRUE)

# only indole with benzyle nitrile has above 
#                 row col
#LNIndole          25  24
#LNBenzyl-nitrile  24  25

x[lower.tri(x, diag = FALSE)] <- NA
x_long <- na.omit(data.frame(melt(x, value.name = "corr", varnames=c('trait_1', 'trait_2'))))

# make heatmap
p <- ggplot(x_long,aes(x=trait_2,y=trait_1,fill= abs(corr))) + 
  geom_tile(show.legend = FALSE) +
  #scale_fill_viridis_c(name = "F Value", option = "C") +
  scale_fill_gradient(low = "white", high = "red") +
  labs(x= NULL, y= NULL) +
  #labs(x= "Population 1", y= "Population 2") +
  geom_text(aes(x=trait_2,y=trait_1, label = round(corr, 2)), color = "black", size = 2) +
  theme_bw(base_size = 10)

# rotate lables, only works here for some reason
p <- p + rotate_x_text(angle = 45) 

resolution <- 100
svglite(paste0("trait_heatmap.svg", sep=""), width = 1080/resolution, height = 1080/resolution)

p

dev.off()