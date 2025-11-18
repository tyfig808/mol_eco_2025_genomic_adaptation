# model to assess how haldane and af are linked
if (!require(tidyverse)) install.packages('tidyverse')
if (!require(vroom)) install.packages('vroom')
if (!require(data.table)) install.packages('data.table')
if (!require(glmmTMB)) install.packages('glmmTMB')
if (!require(DHARMa)) install.packages('DHARMa')
if (!require(fitdistrplus)) install.packages('fitdistrplus')
if (!require(correlation)) install.packages('correlation')
if (!require(car)) install.packages('car')
if (!require(emmeans)) install.packages('emmeans')
if (!require(ggsignif)) install.packages('ggsignif')



out_wd <- "~/mol_eco_2024/2.4_pheno_evo"

setwd(out_wd)

# load haldane
n_counts_af <- fread(sep = "\t", paste0("n_snps_clear_cmh_gwas_haldane_full_rep_soil_sep.tsv"))

# add bee 
n_counts_af$bee<- str_split_i(n_counts_af$treatment, "_", 3)

# look at sums
n_counts_af[, sum(n_snps), by=c("treatment")]

# check dist of variables
descdist(abs(n_counts_af$avg_haldane), discrete = FALSE) # suggests to use beta or gamma
descdist(abs(n_counts_af$n_snps), discrete = TRUE) 

# check if snps affect avg haldane
count_mod <- glmmTMB(abs(avg_haldane) ~ n_snps + (1|treatment) + (1|treatment/trait), family=Gamma(), data = n_counts_af)
S(count_mod)
#            Estimate Std. Error z value Pr(>|z|)    
#(Intercept)  16.4866     1.0679  15.438   <2e-16 ***
#n_snps       -1.0895     0.4398  -2.477   0.0132 *  

Anova(count_mod)
#Response: abs(avg_haldane)
#        Chisq Df Pr(>Chisq)  
#n_snps 6.1359  1    0.01325 *

sim <- simulateResiduals(fittedModel = count_mod, plot = T, n = 100) 