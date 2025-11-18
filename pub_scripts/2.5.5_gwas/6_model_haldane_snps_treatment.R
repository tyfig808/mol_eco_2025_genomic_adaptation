if (!require(tidyverse)) install.packages('tidyverse')
if (!require(vroom)) install.packages('vroom')
if (!require(data.table)) install.packages('data.table')
if (!require(igraph)) install.packages('igraph')
if (!require(bipartite)) install.packages('bipartite')
if (!require(scales)) install.packages('scales')
if (!require(svglite)) install.packages('svglite')
if (!require(geomorph)) install.packages('geomorph')
if (!require(emmeans)) install.packages('emmeans')
if (!require(car)) install.packages('car')
if (!require(glmmTMB)) install.packages('glmmTMB')
if (!require(DHARMa)) install.packages('DHARMa')
if (!require(fitdistrplus)) install.packages('fitdistrplus')

# set wd to place outfile there

out_wd <- "~/mol_eco_2024/2.4_pheno_evo"

setwd(out_wd)

# read in # of snps shared in cmh gwas
n_counts_af <- fread(sep = "\t", paste0("n_snps_clear_cmh_gwas_haldane_full_rep_soil_sep.tsv"))
snps <- n_counts_af[, sum(n_snps), by=c("treatment")]

# read in hald
hald <- n_counts_af[, mean(abs(avg_haldane)), by=c("treatment")]

# 2 Check # snps against haldane  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
descdist(hald$V1, discrete = FALSE)
hald_snps <- glmmTMB(hald$V1 ~ snps$V1, family = beta_family())

Anova(hald_snps)
#         Chisq Df Pr(>Chisq)    
#snps$V1 9.5869  1    0.00196 **
sim_1 <- simulateResiduals(fittedModel = hald_snps, plot = T, n = 100)


# model number snps associated wth phenotype and genomic divergence
setwd("~/mol_eco_2024/2.5.1_fst_pair")
fst <- fread(sep = "\t", paste0("fst_between_g1_and_all_treat.tsv"))


fst_snps <- glmmTMB(fst_full[-1,]$m ~ snps$V1, family = beta_family())
Anova(fst_cmh)

