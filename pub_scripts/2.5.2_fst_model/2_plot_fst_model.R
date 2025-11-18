if (!require(tidyverse)) install.packages('tidyverse')
if (!require(data.table)) install.packages('data.table')
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

# model output
setwd("~/mol_eco_2024/2.5.2_fst_model")

# read in fst and add sim number to move from wide to long, rename cols
mod_fst <- fread(sep = "\t", paste0("fst_rrpp_perm_boxcox_g10_only_trans.tsv"))
mod_fst$sim_num <- 1:nrow(mod_fst)
setnames(mod_fst, c("soil", "herb", "bee", "soil:herb", "soil:bee", "herb:bee"), 
  c("Soil", "Herbivory", "Pollination", "Soil:Herbivory", "Soil:Pollination", "Herbivory:Pollination"),
  skip_absent=TRUE)

# remove the random effects and residuals
col_vec <- c("rep", "chr", "chr:window")
mod_fst <- mod_fst[, (col_vec) := NULL]

# move to long and change name
t_fst = melt(mod_fst, id.vars = c("sim_num"),
                measure.vars = c("Soil", "Herbivory", "Pollination", 
                  "Soil:Herbivory", "Soil:Pollination", "Herbivory:Pollination"))

setnames(t_fst, c("variable", "value"), c("fac_int", "f_value"))

setDT(t_fst)

#t_r <- t_fst[, round(range(f_value, na.rm = T), 2), by = fac_int]
t_m <- t_fst[, round(mean(f_value, na.rm = T), 2), by = fac_int]
t_s <- t_fst[, round(sd(f_value, na.rm = T), 2), by = fac_int]

t_dt <- data.table(cbind(term = levels(t_m$fac_int), mean = t_m$V1, sd = t_s$V1))

# set up plotting window

resolution=144
svglite(paste0("fst_model_g10_rrpp.svg", sep=""), width = 1440/resolution, height = 1080/resolution)
ggplot(t_fst, aes(x = f_value, y = fac_int)) + 
  geom_density_ridges_gradient() +
  coord_cartesian(xlim = c(0,30)) +
  #scale_fill_viridis_c(name = "F Value", option = "C") +
  labs(x = "F Value", y = NULL) + 
  theme_minimal(base_size = 19)

dev.off()


# get percent increase between soil and soil x poll
soil_poll <- t_fst[, mean(f_value), fac_int == "Soil:Pollination"]
soil <- t_fst[, mean(f_value), fac_int == "Soil"]

round((soil_poll[fac_int == TRUE]$V1 - soil[fac_int == TRUE]$V1)/ soil[fac_int == TRUE]$V1, 4)*100
# 168

# get percent increase between herb and herb x poll
herb_poll <- t_fst[, mean(f_value), fac_int == "Herbivory:Pollination"]
#herb <- t_fst[, mean(f_value), fac_int == "Herbivory"] # doesnt work for whatever reason
herb <- mean(t_fst[fac_int == "Herbivory"]$f_value, na.rm = TRUE)

round((herb_poll[fac_int == TRUE]$V1 - herb)/ herb, 4)*100
# 526.64

# get percent increase between soil x herb and soil x herb x poll
soil_herb <- t_fst[, mean(f_value), fac_int == "Soil:Herbivory"]
soil_herb_poll <- t_fst[, mean(f_value), fac_int == "Soil:Herbivory:Pollination"]

round((soil_herb_poll[fac_int == TRUE]$V1 - soil_herb[fac_int == TRUE]$V1)/ soil_herb[fac_int == TRUE]$V1, 4)*100
# 33.17