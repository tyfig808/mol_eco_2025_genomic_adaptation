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

# load fst
# model number snps associated wth phenotype and genomic divergence
setwd("~/mol_eco_2024/2.5.1_fst_pair")
fst <- fread(sep = "\t", paste0("fst_between_g1_and_all_treat.tsv"))

# model if number of selected SNPs associated with phenotypic traits are related to genomic divergence
fst_snps <- glmmTMB(fst_full[-1,]$m ~ snps$V1, family = beta_family())
Anova(fst_snps)

simulateResiduals(fittedModel = fst_snps, plot = T, n = 100)