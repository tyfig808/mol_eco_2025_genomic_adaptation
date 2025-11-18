setwd("~/mol_eco_2024/2.5.3_sel_scans")

if (!require(tidyverse)) install.packages('tidyverse')
if (!require(vroom)) install.packages('vroom')
if (!require(data.table)) install.packages('data.table')
if (!require(poolSeq)) install.packages('poolSeq')

library(compiler)
#if (!require(remotes)) install.packages('remotes')
#remotes::install_github("MartaPelizzola/ACER") # if need from github
library(ACER)

# loop over treatment, can group all together 
treat <- c("LHB", "LHH", "LNHB", "LNHH", "THB", "THH", "TNHB", "TNHH")

# load in effective population size
ne_dt <- fread(sep = ",", paste0("avg_ne.csv"))
#   treatment       V1 bee rounded_ne
#1:       LHB 21.22222   B         22
#2:       LHH 26.77778   H         27
#3:      LNHB 21.77778   B         22
#4:      LNHH 28.55556   H         29
#5:       THB 20.66667   B         21
#6:       THH 27.27778   H         28
#7:      TNHB 22.61111   B         23
#8:      TNHH 28.11111   H         29

## enter the g1 data, subset a then b
freq_a_1 <- vroom("dr_RA_G1.frq", delim = "\t",
           col_names = c("chr", "pos", "nalleles", "nchr_a1", "G1_a1", "G1_a2"), skip = 1) %>% 
          tidyr::unite("ID", c("chr", "pos"), sep="_", remove = FALSE)

freq_b_1 <- vroom("dr_RB_G1.frq", delim = "\t",
           col_names = c("chr", "pos", "nalleles", "nchr_b1", "G1_b1", "G1_b2"), skip = 1) %>% 
          tidyr::unite("ID", c("chr", "pos"), sep="_", remove = FALSE)

# actually faster to find the snps that are present in all then work on the subset
freq_a_1 <- as.data.table(freq_a_1)
freq_b_1 <- as.data.table(freq_b_1)

setkey(freq_a_1, ID)
setkey(freq_b_1, ID)


### create CMH, vroom loads it in better than fread for some reason, slower though
for(i in 1:(length(treat))) {
	subset = treat[i]
	print(paste0("loading files for subset ", subset, "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"))

	## enter the g10 data subset ra_T_H_B then rb_TH_B
	freq_a_10 <- vroom(paste0("dr_RA_", subset, ".frq"), delim = "\t",
	           col_names = c("chr", "pos", "nalleles", "nchr_a10", "G10_a1", "G10_a2"), skip = 1) %>% 
	          tidyr::unite("ID", c("chr", "pos"), sep="_", remove = FALSE)

	freq_b_10 <- vroom(paste0("dr_RB_", subset, ".frq"), delim = "\t",
	           col_names = c("chr", "pos", "nalleles", "nchr_b10", "G10_b1", "G10_b2"), skip = 1) %>% 
	          tidyr::unite("ID", c("chr", "pos"), sep="_", remove = FALSE)
	
	freq_a_10 <- as.data.table(freq_a_10)
	freq_b_10 <- as.data.table(freq_b_10)

	setkey(freq_a_10, ID)
	setkey(freq_b_10, ID)

	# load in g10 prune
	g10_pruned_a <- fread(sep = "\t", header = FALSE, 
			paste0("subset", subset, "_A.prune.in"))
	g10_pruned_b <- fread(sep = "\t", header = FALSE, 
			paste0("subset", subset, "_B.prune.in"))

	# for in file sep by :, here we change to _
	g10_pruned_a$V1 <- sub(':','_', g10_pruned_a$V1)
	setkey(g10_pruned_a, V1)

	g10_pruned_b$V1 <- sub(':','_', g10_pruned_b$V1)
	setkey(g10_pruned_b, V1)

	# merge and take where match in both, add to cmh_list for merging later
	freq_a_10 <- freq_a_10[g10_pruned_a, nomatch = NULL]
	freq_b_10 <- freq_b_10[g10_pruned_b, nomatch = NULL]

	# merge and remove 0
	dtl <- list(freq_a_1, freq_b_1, freq_a_10, freq_b_10)
	B_T_B <- na.omit(Reduce(function(x, y) x[y, on = .(ID)], dtl))

	# allele frequency here is combine so split into two split the G1_a1 col into two
	B_T_B[, c("G1a_base", "G1a_af") := tstrsplit(G1_a1, ":", fixed=TRUE)]
	B_T_B[, c("G1b_base", "G1b_af") := tstrsplit(G1_b1, ":", fixed=TRUE)]
	B_T_B[, c("G10a_base", "G10a_af") := tstrsplit(G10_a1, ":", fixed=TRUE)]
	B_T_B[, c("G10b_base", "G10b_af") := tstrsplit(G10_b1, ":", fixed=TRUE)]

	# make numneric 
	B_T_B$G1a_af <- as.numeric(B_T_B$G1a_af)
	B_T_B$G1b_af <- as.numeric(B_T_B$G1b_af)
	B_T_B$G10a_af <- as.numeric(B_T_B$G10a_af)
	B_T_B$G10b_af <- as.numeric(B_T_B$G10b_af)

	# again na omit with the af and remove where the nchr is 0
	B_T_B <- na.omit(B_T_B, cols=c("G1a_af","G1b_af", "G10a_af", "G10b_af", 
	                    "nchr_a1", "nchr_b1", "nchr_a10", "nchr_b10"))
	B_T_B <- B_T_B[G1a_af != 0 | G1b_af != 0 | G10a_af != 0 | G10b_af != 0]
	B_T_B <- B_T_B[nchr_a1 != 0 | nchr_b1 != 0 | nchr_a10 != 0 | nchr_b10 != 0]

	### perform test on empirical data, format as RepAgen1 RepAgen10 RepBgen1 RepBgen10 
	af_lst <- list(a1 = B_T_B$G1a_af, a10 = B_T_B$G10a_af, b1 = B_T_B$G1b_af, b10 = B_T_B$G10b_af)
	af_mat <- as.matrix(do.call(cbind, af_lst))

	cov_lst <- list(a1 = B_T_B$nchr_a1, a10 = B_T_B$nchr_a10, b1 = B_T_B$nchr_b1, b10 = B_T_B$nchr_b10)
	cov_mat <- as.matrix(do.call(cbind, cov_lst))

	# compute the pvalue with the cmh that takes into account gen and Ne, sig suggests that it is under selection
	# from Thomas science article, the Ne is around 17, rather than 49/2 = 25
	B_T_B$pvalue <- adapted.cmh.test(freq=af_mat, coverage=cov_mat, Ne=rep(ne_dt$rounded_ne[i], 2), gen=c(0,10), repl=1:2, poolSize=NULL)

	# adjust with fdr
	B_T_B$fdr<- p.adjust(B_T_B$pvalue, method = "fdr")
	sub_cmh_fdr <- B_T_B[fdr<0.05] # subset with only sig fdr 

	fwrite(sub_cmh_fdr, sep = "\t", paste0("subset_", subset, "cmh_fdr_sig_updated_ne_pruned_prior.tsv")) # fdr correction and snps sig after
}

