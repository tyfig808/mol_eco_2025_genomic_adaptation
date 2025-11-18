# merge pop sync file with g1 and g1 and for both reps
setwd("~/mol_eco_2024/2.5.3_sel_scans")

if (!require(tidyverse)) install.packages('tidyverse')
if (!require(vroom)) install.packages('vroom')
if (!require(data.table)) install.packages('data.table')


# read in G1 files
g1_a <- fread(sep = "\t", paste0("RA_G1_merged_ref.sync"))
g1_b <- fread(sep = "\t", paste0("RB_G1_merged_ref.sync"))

# loop over treatment
treat <- c("LHB", "LHH", "LNHB", "LNHH", "THB", "THH", "TNHB", "TNHH")

for(i in 1:(length(treat))) {
	s = treat[i]
	print(paste0("loading files for subset ", s, " - - - - - - - - - - - - - - - - - - - - - - - - -"))

	# read in sync file
	sync_dt_a <- fread(sep = "\t", paste0("RA_", s, "_merged_ref.sync"))
	sync_dt_b <- fread(sep = "\t", paste0("RB_", s, "_merged_ref.sync"))

	# combine into treatment sync file
	dtl <- list(g1_a, g1_b, sync_dt_a, sync_dt_b)
	dt <- na.omit(Reduce(function(x, y) x[y, on = .(chr, pos)], dtl))

	# subset and set names 
	dt <- dt[, c("chr", "pos", "ref", "al", "i.al", "i.al.1", "i.al.2")]
	setnames(dt, c("al", "i.al", "i.al.1", "i.al.2"), c("g1_a_al", "g1_b_al", "g10_a_al", "g10_b_al"), skip_absent=TRUE)	

	fwrite(dt, sep = "\t", paste0("treatment_", s, "_g1_g10.sync"))


}





