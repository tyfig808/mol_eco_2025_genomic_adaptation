# create pruned sync file with g1 and g1 and for both reps

if (!require(tidyverse)) install.packages('tidyverse')
if (!require(vroom)) install.packages('vroom')
if (!require(data.table)) install.packages('data.table')

# loop over treatment and create pruned
treat <- c("LHB", "LHH", "LNHB", "LNHH", "THB", "THH", "TNHB", "TNHH")
pruned <- paste(treat, "A", sep = "")


for(i in 1:(length(treat))) {
	s = treat[i]
	print(paste0("loading files for subset ", s, " - - - - - - - - - - - - - - - - - - - - - - - - -"))

	# read in snps sig in cmh with updated ne
	setwd("~/mol_eco_2024/2.5.3_sel_scans")
	dt <- fread(sep = "\t", paste0("treatment_", s, "_g1_g10.sync"))
	setkey(dt, chr, pos)

	# read in list of from plink that are present after ld analyses
	p <- pruned[i]
	clean_snps <- fread(sep = "\t", header = FALSE, paste0("pruned_", p, ".prune.in"))

	# for in file sep by :, here we change to _
	clean_snps$V1 <- sub(':','_', clean_snps$V1)
	
	# create chr and pos to match
	clean_snps$chr = as.character(lapply(strsplit(as.character(clean_snps$V1), split="_"), "[", 1))
	clean_snps$pos = as.numeric(lapply(strsplit(as.character(clean_snps$V1), split="_"), "[", 2))

	setkey(clean_snps, chr, pos)

	# merge and take where match in both, remove the ID col after
	pruned_sync <- dt[clean_snps, nomatch = NULL]
	setkey(pruned_sync, chr, pos)
	pruned_sync[, V1 := NULL]

	# write out
	fwrite(pruned_sync, sep = "\t", col.names = FALSE, # we dont want col names
		paste0("treatment_", s, "_g1_g10_pruned.sync"))

}

