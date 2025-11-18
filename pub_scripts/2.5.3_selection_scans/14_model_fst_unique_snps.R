if (!require(data.table)) install.packages('data.table')
if (!require(GenomicRanges)) install.packages('GenomicRanges')

# load in full cmh file
setwd("~/mol_eco_2024/2.5.3_sel_scans")
clear_total <- fread(sep = "\t", paste0("clear_and_cmh_all_treat_pruned.tsv"))

treat <- c("LHB", "LHH", "LNHB", "LNHH", "THB", "THH", "TNHB", "TNHH")

# make abs col
clear_total[, s_abs := abs(s)]

# create new sig col is sig in cmh and above 0.5 s
clear_total$both_sig <- ifelse(clear_total$s_abs > 0.5 & clear_total$sig == TRUE, 1, 0)

# get count
b_sig <- clear_total[both_sig == 1]


# split into list
sig_list <- split(b_sig, by = 'treatment')

# get the numerical combo list
pair_combos <- combn(c(1:length(sig_list)),2) 

# create list to save results
res_list <- vector("list", length = ncol(pair_combos))

# loop over pairwise combos
for(i in 1:(ncol(pair_combos))) {
  # get names of treat pairs
  treat_1 <- names(sig_list[pair_combos[1,i]])
  treat_2 <- names(sig_list[pair_combos[2,i]])
  print(paste0("run pairwise for treatments: ", treat_1, " - ", treat_2))

  # save dt 
  df1 = sig_list[pair_combos[1,i]][[1]]
  df2 = sig_list[pair_combos[2,i]][[1]]

  # overlap with granges
  gr1 <- with(df1, GRanges(chr, IRanges(start=pos, end=pos, names=ID), SNP=ID))
  gr2 <- with(df2, GRanges(chr, IRanges(start=pos, end=pos, names=ID), SNP=ID))
  ranges <- subsetByOverlaps(gr1, gr2)
  hits <- findOverlaps(gr1, gr2)

  # get subsets 
  sub_1 <- df1[queryHits(hits)]
  sub_2 <- df2[subjectHits(hits)]

  res <- data.table(cbind(ID_1 = sub_1$ID, treatment_1 = treat_1, ID_2 =sub_2$ID, treatment_2 = treat_2))

  # save as res
  res_list[[i]] <- res
  
}

# bind to list, some dont have any matches so fill = true
full_shared <- rbindlist(res_list, fill=TRUE)

# remove the NA and count by pair wise
shared_filt <- na.omit(full_shared)

# get not shared but sig
sig_not_shared <- b_sig[which(!(b_sig$ID %in% shared_filt$ID_1)),]    

n_sig_not_shared <- sig_not_shared[, .N, by = treatment]

# convert to long format to fill easier
shared_long <- data.frame(shared_filt[, .N, by = c("treatment_1", "treatment_2")])

# check if differs with FSt
setwd("~/mol_eco_2024/2.5.1_fst_pair")
fst_full <- fread(sep = "\t", paste0("fst_between_g1_and_all_treat.tsv"))

fst_cmh <- glmmTMB(fst_full[-1,]$m ~ n_sig_not_shared$N, family = beta_family())
Anova(fst_cmh)

#n_sig_not_shared$N 8.2944  1   0.003977 **
simulateResiduals(fittedModel = fst_cmh, plot = T, n = 100) # residual from beta and normal both simular so use beta because it gives smaller p value