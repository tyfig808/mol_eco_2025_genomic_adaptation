# single models to test, fst, n snps, and network size
# run simulations of distance matrix to see if we can say two matrixses are dissimilar
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
file_wd <- "~/mol_eco_2024/2.5.5_gwas/"
out_wd <- "~/mol_eco_2024/2.5.5_gwas/"

# read in files  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
setwd("~/mol_eco_2024/2.5.1_fst_pair")
full_fst <- fread(sep = "\t", paste0("total_rep_sep_fst_full_data.tsv"))

# get average across rep
fst_long <- full_fst[, mean(WEIGHTED_FST), by = c("pop_1", "pop_2", "rep")] 

# take only differences from G1
g1_dt <- fst_long[pop_1 == "G1"]
dt <- g1_dt[, mean(V1), by = pop_2] # take average across reps

# read in network size
setwd(out_wd)
size_dt <- fread(sep = "\t", paste0("network_size_treatments_sigzone_only.tsv"))
size_dt$pop_2 <- paste(size_dt$soil, size_dt$herb, size_dt$bee, sep = "")

# read in shared snps with cmh and clear
setwd("~/mol_eco_2024/2.5.3_sel_scans")
cmh <- fread(sep = "\t", paste0("model_both_clear_cmh_num_snps_soil_sep.tsv"))

#  Check # divergence against netsize  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
descdist(as.numeric(scale(size_dt$size, center = F)), discrete = FALSE)
fst_size <- glmmTMB(1/size_dt$size ~ dt$V1, family = beta_family())
Anova(fst_size)
#       Chisq Df Pr(>Chisq)   
#dt$V1 4.4128  1    0.03567 *

#  Check # n snps against netsize  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
snps_size <- glmmTMB(1/size_dt$size ~ snps$V1, family = beta_family())
Anova(snps_size)
#         Chisq Df Pr(>Chisq)   
#snps$V1 6.8969  1   0.008634 **



