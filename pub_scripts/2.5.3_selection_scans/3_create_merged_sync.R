# merge pop sync file with reference allele
setwd("~/mol_eco_2024/2.5.3_sel_scans")

if (!require(tidyverse)) install.packages('tidyverse')
if (!require(vroom)) install.packages('vroom')
if (!require(data.table)) install.packages('data.table')

# loop over fr files, should be 18
dr_files <- list.files(pattern = "_total.sync", recursive = TRUE, full.names = TRUE)
length(dr_files)==18

# get rep and sub for naming
rep <-  str_split_i(dr_files, "/", 2)
rep <- str_split_i(rep, "_", 1)
sub <-  str_split_i(dr_files, "_", 2)

rep_1 <- str_sub(rep, 2, 2)

for(i in 1:(length(dr_files))) {
	r = rep[i]
	s = sub[i]
	print(paste0("loading files for subset ", r, "_", s, " - - - - - - - - - - - - - - - - - - - - - - - - -"))

	# read in sync file
	sync_dt <- fread(sep = "\t", paste0(r,"_", s, "_total.sync")) 

	# read in ref file, need new chr
	r1 <- rep_1[i]
	ref_dt <- fread(paste0("subset", s,"_", r1, "_chr_pos_ref.tsv")) 
	setnames(ref_dt, c("V1", "V2", "V3"), c("chr", "pos", "ref"), skip_absent=TRUE)	

	# set key and bind
	setkey(sync_dt, chr, pos)
	setkey(ref_dt, chr, pos)
	sync_dt <- sync_dt[ref_dt, nomatch = NULL]

	# remove first ref col and rename
	sync_dt[, ref:= NULL]
	setnames(sync_dt, c("i.ref"), c("ref"), skip_absent=TRUE)	
	setcolorder(sync_dt, c("chr", "pos", "ref", "al"))

	# write out
	fwrite(sync_dt, sep = "\t", paste0(r,"_", s, "_merged_ref.sync"))

}