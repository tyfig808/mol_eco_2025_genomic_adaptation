if (!require(tidyverse)) install.packages('tidyverse')
if (!require(data.table)) install.packages('data.table')
if (!require(ggpubr)) install.packages('ggpubr')
if (!require(viridis)) install.packages('viridis')
if (!require(scales)) install.packages('scales')
if (!require(ggnewscale)) install.packages('ggnewscale')
if (!require(svglite)) install.packages('svglite')
if (!require(igraph)) install.packages('igraph')
if (!require(FNN)) install.packages('FNN')

# read in fst long
setwd("/shares/schiestl.systbot.uzh/variants/data/tfig/ch_1_sub/2.5.1_fst_pair")
full_fst <- fread(sep = "\t", paste0("total_rep_sep_fst_full_data.tsv"))

# get average across rep
fst_long <- full_fst[, mean(WEIGHTED_FST), by = c("pop_1", "pop_2", "rep")] 

# get wide
fst <- dcast(pop_1 + pop_2 ~ rep, value.var = "V1", data = fst_long)

# rename
setnames(fst, c("RA", "RB"), c("RA_fst", "RB_fst"))

# now test with full fst, not just from G1
data2 <- data.frame(pop_1 = fst$pop_2, pop_2 = fst$pop_1, RA_fst = fst$RA_fst) # flip the pop 1 and two 
data3 <- as.data.frame(fst) 
data3<- data3[c("pop_1", "pop_2", "RA_fst")]
df2 <- rbind(data3, data2)
df3 <- as.data.frame.matrix(xtabs(RA_fst ~ ., df2))
rownames(df3) <- colnames(df3)
diag(df3) <- NA

a_mean_fst <- colMeans(df3, na.rm = TRUE)

# convert for RB
data2 <- data.frame(pop_1 = fst$pop_2, pop_2 = fst$pop_1, RB_fst = fst$RB_fst) # flip the pop 1 and two 
data3 <- as.data.frame(fst) 
data3<- data3[c("pop_1", "pop_2", "RB_fst")]
df2 <- rbind(data3, data2)
df3 <- as.data.frame.matrix(xtabs(RB_fst ~ ., df2))
rownames(df3) <- colnames(df3)
diag(df3) <- NA

b_mean_fst <- colMeans(df3, na.rm = TRUE)

# get mean
mean_fst_g1 <- data.table(cbind.data.frame(treat = names(b_mean_fst), a_fst = a_mean_fst, b_fst = b_mean_fst))
mean_fst_g1[, m := rowMeans(.SD), by = treat]

fwrite(mean_fst_g1, sep = "\t", paste0("fst_between_g1_and_all_treat.tsv"))

# add bee 
mean_fst_g1$bee <- str_sub(mean_fst_g1$treat, -1)

# run model
bee_fst <- glmmTMB(m~bee, mean_fst_g1[-1])
Anova(bee_fst)
#     Chisq Df Pr(>Chisq)    
#bee 55.512  1  9.287e-14 ***

mean_fst_g1[, mean(m), by = bee]
mean_fst_g1[, sd(m), by = "bee"]

#2:   B 0.10253126
#3:   H 0.07877896

#2:   B 0.005519437
#3:   H 0.004871764
