# create sync file for clear analyses
setwd("~/mol_eco_2024/2.5.3_sel_scans")

if (!require(tidyverse)) install.packages('tidyverse')
if (!require(vroom)) install.packages('vroom')
if (!require(data.table)) install.packages('data.table')

# get freq files from PLINK2
# loop over frq files, should be 18
dr_files <- list.files(pattern = "dr.*fr", recursive = TRUE, full.names = TRUE)
length(dr_files)==18

# get rep and sub for naming
rep <- str_split_i(dr_files, "_", 2)
sub <-  str_split_i(dr_files, "_", 3)
sub <-  str_split_i(sub, "[.]", 1)

for(i in 1:(length(dr_files))) {
	r = rep[i]
	s = sub[i]
	print(paste0("loading files for subset ", r, "_", s, " - - - - - - - - - - - - - - - - - - - - - - - - -"))

	## enter the g1 data, subset a then b
	freq_a_1 <- vroom(dr_files[i], delim = "\t",
	           col_names = c("chr", "pos", "nalleles", "nchr_a1", "G1_a1", "G1_a2"), skip = 1) %>% 
	          tidyr::unite("ID", c("chr", "pos"), sep="_", remove = FALSE)

  	print(paste0("done loading files for subset ", r, "_", s, " - - - - - - - - - - - - - - - - - - - - - - - - -"))

	# conver to data table and seperate base from freq
	freq_a_1 <- data.table(freq_a_1)
	freq_a_1[, c("G1a_base", "G1a_af") := tstrsplit(G1_a1, ":", fixed=TRUE)]
	freq_a_1[, c("G1a2_base", "G1a2_af") := tstrsplit(G1_a2, ":", fixed=TRUE)]

	# make numeric
	freq_a_1$G1a_af <- as.numeric(freq_a_1$G1a_af)
	freq_a_1$G1a2_af <- as.numeric(freq_a_1$G1a2_af)

	# remove na
	freq_a_1 <- na.omit(freq_a_1, cols=c("G1a_af","G1a2_af", "nchr_a1"))

	# get counts per allele and round to whole number
	freq_a_1[, a1_count := round(nchr_a1 * G1a_af, 0)]
	freq_a_1[, a2_count := round(nchr_a1 * G1a2_af, 0)]

	# check to see if they sum up, we want to see 0
	sum(which(freq_a_1$a1_count + freq_a_1$a2_count != freq_a_1$nchr_a1))

	# we want to make this a sync file,
	#2R  2302    N   0:7:0:0:0:0 0:7:0:0:0:0
	#2R  2303    N   0:8:0:0:0:0 0:8:0:0:0:0
	#2R  2304    N   0:0:9:0:0:0 0:0:9:0:0:0
	#2R  2305    N   1:0:9:0:0:0 0:0:9:1:0:0

	#col1: reference contig
	#col2: position within the refernce contig
	#col3: reference character
	#col4: allele frequencies of population number 1
	#col5: allele frequencies of population number 2
	#coln: allele frequencies of population number n
	# like so: in the format A:T:C:G:N:del, 

	# create data table of same length
	v <- vector(length = nrow(freq_a_1), "numeric")
	dt <- data.table("A" = v, "T" = v, "C" = v, "G" = v, "N" = v, "del" = v)  # Create data.table

	# check, if base is the base specified, put the counts, otherwise 0
	dt$A <- ifelse(freq_a_1$G1a_base == "A", freq_a_1$a1_count, 0)
	dt$T <- ifelse(freq_a_1$G1a_base == "T", freq_a_1$a1_count, 0)
	dt$C <- ifelse(freq_a_1$G1a_base == "C", freq_a_1$a1_count, 0)
	dt$G <- ifelse(freq_a_1$G1a_base == "G", freq_a_1$a1_count, 0)

	# add the second allele, leave as is otherwise
	dt$A <- ifelse(freq_a_1$G1a2_base == "A", freq_a_1$a2_count, dt$A)
	dt$T <- ifelse(freq_a_1$G1a2_base == "T", freq_a_1$a2_count, dt$T)
	dt$C <- ifelse(freq_a_1$G1a2_base == "C", freq_a_1$a2_count, dt$C)
	dt$G <- ifelse(freq_a_1$G1a2_base == "G", freq_a_1$a2_count, dt$G)

	# check to see if only two non zero rows, should be 0 
	which(rowSums(dt != 0)!= 2)

	# create paste together dt to get the count column we want
	al <- dt[,  do.call(paste, c(.SD, sep = ":"))]

	# save as sync file and write out
	sync_dt <- cbind(chr = freq_a_1$chr, pos = freq_a_1$pos, ref = v, al = al)

	fwrite(sync_dt, sep = "\t", paste0(r,"_", s, "_total.sync")) 


}


