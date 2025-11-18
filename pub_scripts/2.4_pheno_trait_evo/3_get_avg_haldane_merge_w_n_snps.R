# add haldane to shared af cmh file
setwd("~/mol_eco_2024/2.4_pheno_evo")

if (!require(tidyverse)) install.packages('tidyverse')
if (!require(data.table)) install.packages('data.table')
if (!require(stringr)) install.packages('stringr')
if (!require(stringi)) install.packages('stringi')

# load in 
full_shared <- fread(sep = "\t", paste0("shared_cmh_af_all_treatments_traits.tsv")) 
full_hald <- fread(sep = "\t", paste0("haldane_full_rep_soil_sep.tsv"))

# get mean haldane for treatment
mean_haldane <- full_hald[, mean(haldane), by = c("treatment", "trait")]
setnames(mean_haldane, c("V1"), c("avg_haldane"), skip_absent=TRUE)	

# change treatment names in full_shared, add underscores to these positions
stri_sub(full_shared$treatment, 2, 1) <- "_"
stri_sub(full_shared$treatment, -1, -2) <- "_"

# merge
setkey(full_shared, treatment, trait)
setkey(mean_haldane, treatment, trait)
af_hald <- mean_haldane[full_shared]

# write out
fwrite(af_hald, sep = "\t", paste0("shared_af_haldane_full_rep_soil_sep.tsv"))

# get number of traits per treatment, trait, merge with hald
n_counts_af <- af_hald[, sum(!is.na(trans_af_diff)), by=c("treatment", "trait")]

setkey(n_counts_af, treatment, trait)
setkey(mean_haldane, treatment, trait)
n_counts_af <- mean_haldane[n_counts_af]
setnames(n_counts_af, c("V1"), c("n_snps"), skip_absent=TRUE)	

fwrite(n_counts_af, sep = "\t", paste0("n_snps_af_haldane_full_rep_soil_sep.tsv"))