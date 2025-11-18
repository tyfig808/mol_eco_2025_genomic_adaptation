if (!require(data.table)) install.packages('data.table')
if (!require(igraph)) install.packages('igraph')
if (!require(bipartite)) install.packages('bipartite')
if (!require(emmeans)) install.packages('emmeans')
if (!require(ggsignif)) install.packages('ggsignif')
if (!require(car)) install.packages('car')
if (!require(effects)) install.packages('effects')
if (!require(ggridges)) install.packages('ggridges')
if (!require(glmmTMB)) install.packages('glmmTMB')
if (!require(DHARMa)) install.packages('DHARMa')
if (!require(UpSetR)) install.packages('UpSetR')
library(stringr)

# test to see if density is related to the number of between trait groups

# set wd to place outfile there
out_wd <- "~/mol_eco_2024/2.5.5_gwas/"

# load trait matrix full file
setwd(out_wd)
t_mat <- fread(sep = "\t", paste0("shared_snps_between_traits_sigzone_only.tsv"))

# load in size
size_dt <- fread(sep = "\t", paste0("network_size_treatments_sigzone_only.tsv"))

# add merged trait column
t_mat[, pair_trait := paste(trait_1, trait_2, sep = "_")]

# get number of shared traits between classes
t <- t_mat[, sum(snps_shared), by = c("treat", "shared_class")]

# test to see how shared traits affect size -------------------------------------------------------------
m_v_dt <- t[shared_class == "morphology_volatile"]

# model to see if density a function of shared traits between classes
m_v_mod <- glmmTMB(1/size_dt$size ~ m_v_dt$V1, family = beta_family())
Anova(m_v_mod)

#           Chisq Df Pr(>Chisq)    
#m_v_dt$V1 25.055  1  5.572e-07 ***