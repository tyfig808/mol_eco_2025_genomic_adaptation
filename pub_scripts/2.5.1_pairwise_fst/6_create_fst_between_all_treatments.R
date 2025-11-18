if (!require(tidyverse)) install.packages('tidyverse')
if (!require(vroom)) install.packages('vroom')
if (!require(data.table)) install.packages('data.table')
if (!require(car)) install.packages('car')
if (!require(glmmTMB)) install.packages('glmmTMB')
if (!require(DHARMa)) install.packages('DHARMa')
if (!require(fitdistrplus)) install.packages('fitdistrplus')

setwd("/shares/schiestl.systbot.uzh/variants/data/tfig/ch_1_sub/2.5.1_fst_pair")
fst_full <- fread(sep = "\t", paste0("total_rep_sep_fst_full_data.tsv"))

fst <- fst_full[, mean(WEIGHTED_FST), by = c("pop_1", "pop_2", "rep")]

fst_a <- fst[rep == "RA"]
setnames(fst_a, "V1", "RA_fst")

data2 <- data.frame(pop_1 = fst_a$pop_2, pop_2 = fst_a$pop_1, RA_fst = fst_a$RA_fst) # flip the pop 1 and two 
data3 <- as.data.frame(fst_a) 
data3<- data3[c("pop_1", "pop_2", "RA_fst")]
df2 <- rbind(data3, data2)
df3 <- as.data.frame.matrix(xtabs(RA_fst ~ ., df2))
rownames(df3) <- colnames(df3)
diag(df3) <- NA

a_mean_fst <- colMeans(df3, na.rm = TRUE)

# convert for RB
fst_b <- fst[rep == "RB"]
setnames(fst_b, "V1", "RB_fst")


data2 <- data.frame(pop_1 = fst_b$pop_2, pop_2 = fst_b$pop_1, RB_fst = fst_b$RB_fst) # flip the pop 1 and two 
data3 <- as.data.frame(fst_b) 
data3<- data3[c("pop_1", "pop_2", "RB_fst")]
df2 <- rbind(data3, data2)
df3 <- as.data.frame.matrix(xtabs(RB_fst ~ ., df2))
rownames(df3) <- colnames(df3)
diag(df3) <- NA

b_mean_fst <- colMeans(df3, na.rm = TRUE)

mean_fst_g1 <- data.table(cbind.data.frame(treat = names(b_mean_fst), a_fst = a_mean_fst, b_fst = b_mean_fst))
mean_fst_g1[, m := rowMeans(.SD), by = treat]

# write out
fwrite(mean_fst_g1, sep = "\t", paste0("fst_between_g1_and_all_treat.tsv"))