if (!require(tidyverse)) install.packages('tidyverse')
if (!require(vroom)) install.packages('vroom')
if (!require(data.table)) install.packages('data.table')
if (!require(emmeans)) install.packages('emmeans')
if (!require(car)) install.packages('car')
if (!require(glmmTMB)) install.packages('glmmTMB')
if (!require(DHARMa)) install.packages('DHARMa')
if (!require(fitdistrplus)) install.packages('fitdistrplus')


# read in files  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
setwd("~/mol_eco_2024/2.5.1_fst_pair")
fst_full <- fread(sep = "\t", paste0("fst_between_g1_and_all_treat.tsv"))

setwd("~/mol_eco_2024/2.5.3_sel_scans")
cmh <- fread(sep = "\t", paste0("model_both_clear_cmh_num_snps_soil_sep.tsv"))

# 1b check snps against total fst, no pheno - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
descdist(fst_full$m, discrete = FALSE) # although says to use beta, residuals better without

fst_cmh <- glmmTMB(fst_full[-1,]$m ~ cmh$num_snps)
Anova(fst_cmh)
#		  Chisq Df Pr(>Chisq)     
#cmh$num_snps 19.264  1   1.139e-05 ***
sim_1 <- simulateResiduals(fittedModel = fst_cmh, plot = T, n = 100)