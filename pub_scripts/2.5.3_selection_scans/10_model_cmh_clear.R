# run model on clear s, use abs and maybe need to transform
if (!require(data.table)) install.packages('data.table')
if (!require(tidyverse)) install.packages('tidyverse')
if (!require(glmmTMB)) install.packages('glmmTMB')
if (!require(DHARMa)) install.packages('DHARMa')
if (!require(fitdistrplus)) install.packages('fitdistrplus')
if (!require(car)) install.packages('car')
if (!require(emmeans)) install.packages('emmeans')
if (!require(gammit)) install.packages('gammit')


library(mgcv) # use this with bam
require(parallel) # this may speed up since we have big boys

x <- detectCores()
nc <- x/(8/7)  ## use 80% of available cores, set for example portability
ctrl <- list(nthreads=nc)

if (detectCores()>1) { ## no point otherwise
  cl <- makeCluster(nc) 
  ## could also use makeForkCluster, but read warnings first!
} else cl <- NULL

# load in full cmh file
setwd("~/mol_eco_2024/2.5.3_sel_scans")
clear_total <- fread(sep = "\t", paste0("clear_and_cmh_all_treat_pruned.tsv"))

# make treatment cols
clear_total$soil <- str_sub(clear_total$treatment, 1, 1)
clear_total$herb <- str_sub(clear_total$treatment, 2, -2)
clear_total$bee <- str_sub(clear_total$treatment, -1)

# make abs col
clear_total[, s_abs := abs(s)]

# create new sig col is sig in cmh and above 0.5 s
clear_total$both_sig <- ifelse(clear_total$s_abs > 0.5 & clear_total$sig == TRUE, 1, 0)

# get count
b_sig <- clear_total[both_sig == 1]
treat_sig <- b_sig[, .N, by = c("soil", "herb", "bee")]
setkey(treat_sig, "soil", "herb", "bee")
#   soil herb bee   N
#1:    L    H   B 161
#2:    L    H   H 253
#3:    L   NH   B 729
#4:    L   NH   H  98
#5:    T    H   B 568
#6:    T    H   H  80
#7:    T   NH   B 461
#8:    T   NH   H  74


# total model
clear_total$chr <- as.factor(clear_total$chr)

m1 <- bam(
	both_sig ~  soil * herb * bee + 
 		s(chr, bs ='re') + s(chr, pos, bs ='re'), # random effects
  	data = clear_total, family = "binomial",
  	discrete = TRUE, chunk.size=500000, control=ctrl)

anova(m1)
#              df  Chi.sq p-value
#soil           1 212.899  <2e-16
#herb           1 291.465  <2e-16
#bee            1   5.777  0.0162
#soil:herb      1 297.390  <2e-16
#soil:bee       1 258.104  <2e-16
#herb:bee       1 291.506  <2e-16
#soil:herb:bee  1 139.060  <2e-16

#Approximate significance of smooth terms:
#              edf Ref.df Chi.sq p-value
#s(chr)      8.590 10.000 431864       1
#s(pos,chr)  9.338 10.000  24811       1


b_m <- clear_total[, mean(both_sig), by = c("soil", "herb", "bee")]
b_se <- clear_total[, sd(both_sig)/.N, by = c("soil", "herb", "bee")]
clear_total[, .N, by = c("soil", "herb", "bee")]

dt <- cbind(b_m, se = b_se$V1)
setnames(dt, "V1",  "mean")
setkey(dt, "soil", "herb", "bee")

tot_em <- emmeans(m1,  ~ soil * herb * bee , type = "response")
comp_dt <- data.table(as.data.frame(multcomp::cld(tot_em))) 
comp_dt$group <- as.numeric(comp_dt$.group)
# flip the order because for whatever reason it is saying that the top one is rank 3, flip so rank 1 instead
comp_dt$group <- match(comp_dt$group, sort(unique(comp_dt$group), decreasing = TRUE)) 
setkey(comp_dt, "soil", "herb", "bee")


# run for tuff and lime sep
tuff_cmh <- clear_total[soil == "T"]
t1 <- bam(
	both_sig ~   herb * bee + 
 		s(chr, bs ='re') + s(chr, pos, bs ='re'), # random effects
  	data = tuff_cmh, family = "binomial",
  	discrete = TRUE, chunk.size=500000, control=ctrl)

anova(t1)
#         df  Chi.sq  p-value
#herb      1  34.603 4.04e-09
#bee       1 363.307  < 2e-16
#herb:bee  1   0.932    0.334

#Approximate significance of smooth terms:
#              edf Ref.df Chi.sq p-value
#s(chr)      8.481 10.000 235807       1
#s(pos,chr)  9.327 10.000  14197       1

summary(t1)
#(Intercept) -5.97658    0.29139 -20.511  < 2e-16 ***
#herbNH      -0.36972    0.06285  -5.882 4.04e-09 ***
#beeH        -2.27788    0.11951 -19.061  < 2e-16 ***
#herbNH:beeH  0.16712    0.17314   0.965    0.334    


#
lime_cmh <- clear_total[soil == "L"]
l1 <- bam(
	both_sig ~   herb * bee + 
 		s(chr, bs ='re') + s(chr, pos, bs ='re'), # random effects
  	data = lime_cmh, family = "binomial",
  	discrete = TRUE, chunk.size=500000, control=ctrl)

anova(l1)
#         df  Chi.sq p-value
#herb      1 288.217  <2e-16
#bee       1   5.794  0.0161
#herb:bee  1 290.256  <2e-16

#Approximate significance of smooth terms:
#              edf Ref.df Chi.sq p-value
#s(chr)      8.079 10.000  95367       1
#s(pos,chr)  7.164 10.000   3317       1

summary(l1)
#            Estimate Std. Error z value Pr(>|z|)    
#(Intercept) -7.08401    0.18972 -37.339   <2e-16 ***
#herbNH       1.48004    0.08718  16.977   <2e-16 ***
#beeH         0.24284    0.10088   2.407   0.0161 *  
#herbNH:beeH -2.51355    0.14754 -17.037   <2e-16 ***

# create full soil specific data set
t_em <- emmeans(t1,  ~ herb * bee , type = "response")
t_dt <- data.table(as.data.frame(multcomp::cld(t_em))) 
t_dt$group <- as.numeric(t_dt$.group)
# flip the order because for whatever reason it is saying that the top one is rank 3, flip so rank 1 instead
t_dt$group <- match(t_dt$group, sort(unique(t_dt$group), decreasing = TRUE)) 
t_dt$soil <- "T"

# for limestone now
l_em <- emmeans(l1,  ~ herb * bee , type = "response")
l_dt <- data.table(as.data.frame(multcomp::cld(l_em))) 
l_dt$group <- as.numeric(l_dt$.group)
# flip the order because for whatever reason it is saying that the top one is rank 3, flip so rank 1 instead
l_dt$group <- match(l_dt$group, sort(unique(l_dt$group), decreasing = TRUE)) 
l_dt$soil <- "L"

soil_dt <- rbindlist(list(t_dt, l_dt))
setkey(soil_dt, "soil", "herb", "bee")

dt <- cbind(dt, full_group = comp_dt$group, soil_group = soil_dt$group, num_snps = treat_sig$N)


# get ci from the mean
b_p <- clear_total[, mean(both_sig), by = c("soil", "herb", "bee")]
b_n <- clear_total[, .N, by = c("soil", "herb", "bee")]

b_p$margin <- qnorm(0.975)*sqrt(b_p$V1*(1-b_p$V1)/b_n$N)

b_p[, lcl :=  V1 - margin]
b_p[, ucl :=  V1 + margin]

setkey(b_p, "soil", "herb", "bee")

dt <- cbind(dt, lcl = b_p$lcl, ucl = b_p$ucl)

fwrite(dt, sep = "\t", paste0("model_both_clear_cmh_num_snps_soil_sep.tsv"))