# test overlap with plus minus some kb
if (!require(data.table)) install.packages('data.table')
if (!require(GenomicRanges)) install.packages('GenomicRanges')

# load in full cmh file
setwd("~/mol_eco_2024/2.5.3_sel_scans")
clear_total <- fread(sep = "\t", paste0("clear_and_cmh_all_treat_pruned.tsv"))


# make abs col
clear_total[, s_abs := abs(s)]

# create new sig col is sig in cmh and above 0.5 s
clear_total$both_sig <- ifelse(clear_total$s_abs > 0.5 & clear_total$sig == TRUE, 1, 0)

# get count
b_sig <- clear_total[both_sig == 1]

# create plus and minus with 0 extra, this can be used to check if more snps if larger KB increase, did not increase much with .5 KB extra
kb = 1 * 0

b_sig[, end := pos+kb]
b_sig[, start := pos-kb]

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
  gr1 <- with(df1, GRanges(chr, IRanges(start=start, end=end, names=ID), SNP=ID))
  gr2 <- with(df2, GRanges(chr, IRanges(start=start, end=end, names=ID), SNP=ID))
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

# convert to long format to fill easier
shared_long <- data.frame(shared_filt[, .N, by = c("treatment_1", "treatment_2")])

# create matrix to fill, we need this because not all treatments have all the treats
treat <- c("LHB", "LHH", "LNHB", "LNHH", "THB", "THH", "TNHB", "TNHH")
mat_1 <- matrix(0, nrow = length(treat), ncol = length(treat))
row.names(mat_1) <- treat
colnames(mat_1) <- treat

# fill with long treat format
mat_1[as.matrix(shared_long[c("treatment_1", "treatment_2")])] <- as.double(shared_long$N)

# transpose, if want to remove the diag also
mat_1 <-t(mat_1)

# write out
setwd("~/mol_eco_2024/2.5.3_sel_scans")
write.table(mat_1, file = "shared_snp_exact_pos_mat.csv", sep = ",")