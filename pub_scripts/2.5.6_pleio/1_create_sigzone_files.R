# extract the peak of the sigzones for each trait within treatment, and create a zone of 1kb before and after
if (!require(tidyverse)) install.packages('tidyverse')
if (!require(vroom)) install.packages('vroom')
if (!require(data.table)) install.packages('data.table')

# create for later use	
file_wd <- "~/mol_eco_2024/2.5.5_gwas/"
out_wd <- "~/mol_eco_2024/2.5.5_gwas/"

# this is the kb range used
r <- 1000
kb <- r/1000

### treatments to loop over
treat <- c("LHB", "LHH", "LNHB", "LNHH", "THB", "THH", "TNHB", "TNHH")

for(j in 1:length(treat)) {

	# set wd for treatment
	subset = treat[j]
	print(paste0("loading files for subset ", subset))
	out_dir <- paste(out_wd, subset, "/", sep = "")
	setwd(out_dir)
	#getwd()

	### read in sigzone 
	sz_files <- list.files(pattern = "_emmax_sigzones.tsv", recursive = TRUE, full.names = TRUE)
	t_names <- str_split_i(sz_files, "_", 2)

	for(i in 1:length(sz_files)) {
	#i <- 1	
		trait <- t_names[i]
		print(paste0("loading files for subset ", subset, " and for trait: ", trait))

		# read in files for trait
		sigZones <- fread(sz_files[i])
		
		# add plus minus 1kb at peak and add sigzone number 
		sigZones[, minus_kb := end - r]
		sigZones[, plus_kb := end + r]
		sigZones[, zone_id := paste(CHROM, end, sep = "_")]

		sigZones[, treatment := subset]
		sigZones[, trait := trait]

		# order the sigzone 
		keycol <-c("CHROM","beg")
		sigZones<- setorderv(sigZones, keycol)

	    # write sizone with kb range
	    fwrite(sigZones, sep = "\t", paste0(subset, "_", trait, "_emmax_only_sigzones_", kb, "kb_range.tsv")) #cmh local score with sig zones here 
    }

} 