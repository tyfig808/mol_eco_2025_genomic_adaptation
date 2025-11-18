# create total clear tsv
setwd("~/mol_eco_2024/2.5.3_sel_scans")

if (!require(tidyverse)) install.packages('tidyverse')
if (!require(vroom)) install.packages('vroom')
if (!require(data.table)) install.packages('data.table')

# loop over treatment
treat <- c("LHB", "LHH", "LNHB", "LNHH", "THB", "THH", "TNHB", "TNHH")

clear_list <- vector("list", length = length(treat))

for(i in 1:(length(treat))) {
	subset = treat[i]
	print(paste0("loading files for subset ", subset, " - - - - - - - - - - - - - - - - - - - - - - - - -"))

	# read in sync file
	clear_dt <- fread(sep = "\t", col.names = c("chr", "pos", "alt", "s", "null"), paste0("clear_", subset, ".tsv"))
	
	# comes with three additonal rows so we can remove those
	clear_dt <- clear_dt[4:nrow(clear_dt),] 

	# add treatment col and bind together
	clear_dt[, treatment:= subset]
	
	# add to list
	clear_list[[i]] <- clear_dt

}

clear_total <- rbindlist(clear_list)

fwrite(clear_total, sep = "\t", paste0("clear_all_treat_pruned.tsv"))