if (!require(tidyverse)) install.packages('tidyverse')
if (!require(vroom)) install.packages('vroom')
if (!require(data.table)) install.packages('data.table')
if (!require(car)) install.packages('car')
if (!require(glmmTMB)) install.packages('glmmTMB')
if (!require(DHARMa)) install.packages('DHARMa')
if (!require(fitdistrplus)) install.packages('fitdistrplus')

setwd("/shares/schiestl.systbot.uzh/variants/data/tfig/ch_1_sub/2.5.1_fst_pair")
full_fst <- fread(sep = "\t", paste0("total_rep_sep_fst_full_data.tsv"))

# get average across rep
fst_long <- full_fst[, mean(WEIGHTED_FST), by = c("pop_1", "pop_2", "rep")] 


# take only treatment plants
fst_b <- fst_long[pop_1 != "G1"]

# get the bee
fst_b$bee_1 <- str_sub(fst_b$pop_1, -1)
fst_b$bee_2 <- str_sub(fst_b$pop_2, -1)

fst_b[, b_b := paste(bee_1, bee_2, sep = "_")]
fst_b$b_b <- ifelse(fst_b$b_b == "H_B", "B_H", fst_b$b_b)

b_dt <- fst_b[b_b == "B_B" | b_b == "B_H"]

# check whether b_b larger than h_h
descdist(b_dt$V1, discrete = FALSE)
b_mod <- glmmTMB(V1 ~ b_b + (1|rep), data = b_dt, family = beta_family())

Anova(b_mod)
#     Chisq Df Pr(>Chisq)    
#b_b 93.19  1  < 2.2e-16 ***
simulateResiduals(fittedModel = b_mod, plot = T, n = 100)

# get the stats
b_m <- b_dt[, mean(V1), by = b_b]
b_sd <- b_dt[, sd(V1), by = b_b]

#> b_m
#   b_b         V1
#1: B_H 0.09596604
#2: B_B 0.12359821
#> b_sd
#   b_b          V1
#1: B_H 0.009987473
#2: B_B 0.012566427

