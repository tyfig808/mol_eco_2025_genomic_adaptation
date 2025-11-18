# create total clear tsv
setwd("~/mol_eco_2024/2.5.3_sel_scans")

if (!require(tidyverse)) install.packages('tidyverse')
if (!require(vroom)) install.packages('vroom')
if (!require(data.table)) install.packages('data.table')

# loop over treatment
treat <- c("LHB", "LHH", "LNHB", "LNHH", "THB", "THH", "TNHB", "TNHH")

# cols to take for af
myVector <- c("ID", "chr", "pos", "fdr", "s")

clear_list <- vector("list", length = length(treat))

for(i in 1:(length(treat))) {
	subset = treat[i]
	print(paste0("loading files for subset ", subset, " - - - - - - - - - - - - - - - - - - - - - - - - -"))

	# read in total cmh with updated ne
	cmh <- fread(sep = "\t", paste0("subset_", subset, "cmh_total_updated_ne.tsv")) 
	setkey(cmh, ID)

	# read in clear which is pruned already, # comes with three additonal rows so we can remove those
	clear_dt <- fread(sep = "\t", col.names = c("chr", "pos", "alt", "s", "null"), paste0("clear_", subset, ".tsv"))
	clear_dt <- clear_dt[4:nrow(clear_dt),] 

	# add id col and bind together
	clear_dt[, ID:= paste(chr, pos, sep = "_")]
	setkey(clear_dt, ID)
	clear_cmh <- cmh[clear_dt, nomatch = NULL]

	# take certain cols
	clear_cmh <- clear_cmh[, ..myVector]
	
	# add sig col and treatment col
	clear_cmh[, sig:= ifelse(fdr < 0.05, TRUE, FALSE)]
	clear_cmh[, treatment:= subset]

	# write ind and also add to list
	fwrite(clear_cmh, sep = "\t", paste0("subset_", subset, "_clear_and_cmh_pruned.tsv"))
	clear_list[[i]] <- clear_cmh

}


clear_total <- rbindlist(clear_list)

fwrite(clear_total, sep = "\t", paste0("clear_and_cmh_all_treat_pruned.tsv"))


# create bed file of sig snps
clear_total <- fread(sep = "\t", paste0("clear_and_cmh_all_treat_pruned.tsv"))

# create total bed
total_bed <- clear_total[, c("chr", "pos")]
setkey(total_bed, "chr", "pos")

setnames(total_bed, c("pos"), c("chromStart"), skip_absent=TRUE)	
total_bed$chromEnd <- total_bed$chromStart

fwrite(total_bed, sep = "\t", paste0("clear_cmh_all_treat_pruned.bed"))