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
file_wd <- "~/mol_eco_2024/2.5.5_gwas/"
out_wd <- "~/mol_eco_2024/2.5.5_gwas/"

# load trait matrix full file
setwd(out_wd)
t_mat <- fread(sep = "\t", paste0("shared_snps_between_traits_sigzone_only.tsv"))

# load in size
size_dt <- fread(sep = "\t", paste0("network_size_treatments_sigzone_only.tsv"))

# add merged trait column
t_mat[, pair_trait := paste(trait_1, trait_2, sep = "_")]

m_v_mat <- t_mat[shared_class == "morphology_volatile"]

m_v_mat$bee <- str_sub(m_v_mat$treat, -1)

m_v_mod_bee  <- glmmTMB(snps_shared ~ bee + (1|pair_trait), data = m_v_mat,
		family = "truncated_nbinom2", ziformula= ~ (1|pair_trait))

Anova(m_v_mod_bee)
#     Chisq Df Pr(>Chisq)    
#bee 11.818  1  0.0005866 ***

S(m_v_mod_bee)
