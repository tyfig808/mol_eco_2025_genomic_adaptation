# using sigzones with plus kb and minus get overlap between traits and plot only those effects of snps that occur in both traits
if (!require(tidyverse)) install.packages('tidyverse')
if (!require(vroom)) install.packages('vroom')
if (!require(data.table)) install.packages('data.table')
if (!require(GenomicRanges)) install.packages('GenomicRanges')

#install.packages("remotes")
#remotes::install_github("jokergoo/ComplexHeatmap")
if (!require(ComplexHeatmap)) install.packages('ComplexHeatmap')

# set for later
file_wd <- "~/mol_eco_2024/2.5.5_gwas/"
out_wd <- "~/mol_eco_2024/2.5.5_gwas/"

# take 1 kb around peak
r <- 1000
kb <- r/1000

# load treatments for loop over
treat <- c("LHB", "LHH", "LNHB", "LNHH", "THB", "THH", "TNHB", "TNHH")

# loop over treatment, about 1-2 minutes per treatment
for(j in 1:length(treat)) {

	### treatments to array over
	subset = treat[j]
	print(paste0("loading files for subset ", subset, "- - - - - - - - - - - - - - - - - - - - - - - - -"))

	# set wd to place outfile there
	out_dir <- paste(out_wd, subset, "/", sep = "")
	setwd(out_dir)

	### get sigzone files with beg and end and full local score with distances
	kb_files <- list.files(pattern = "_emmax_only_sigzones_", recursive = TRUE, full.names = TRUE)

	# get trait names and kb size
	t_names <- str_split_i(kb_files, "_", 2)

	# create list for upset plot
	id_list <- vector("list", length = length(t_names))
	trait_list <- vector("list", length = length(t_names))

	# loop over id list, extract those that are within both 
	n <- choose(length(t_names),2) # number of combinations


	name <- vector(length = n)
	res_list <- vector("list", length = n)

	for(i in 1:length(kb_files)) {
		#i = 1
		trait = t_names[i]
		print(paste0("loading files for subset ", subset, " and for trait: ", trait))

		# load in sigzones and create id col
		x <- fread(sep = "\t", paste0(subset, "_", trait, "_emmax_only_sigzones_", kb, "kb_range.tsv"))
		id_list[[i]] <- x$zone_id
	}

	counter <- 1
	for(i in 1:(length(id_list)-1)) {

		# load first trait in pairwise and read in sigzone
		trait = t_names[i]
		sigZones_1 <- fread(sep = "\t", paste0(subset, "_", trait, "_emmax_only_sigzones_", kb, "kb_range.tsv"))

		# create first genomic range object
		gr1 <- GRanges(sigZones_1$CHROM, IRanges(start=sigZones_1$minus_kb, end=sigZones_1$plus_kb), names = sigZones_1$zone_id)
		

		# second loop 
	  	for(k in (i+1):length(id_list)) {

  			# load second trait in pairwise and load second sigzone
			trait2 = t_names[k]
			print(paste0("loading files for subset ", subset, " pariwse between ", trait, " and for ", trait2)) 
			sigZones_2 <- fread(sep = "\t", paste0(subset, "_", trait2, "_emmax_only_sigzones_", kb, "kb_range.tsv"))

			# create second genomic range object
			gr2 <- GRanges(sigZones_2$CHROM, IRanges(start=sigZones_2$minus_kb, end=sigZones_2$plus_kb), names = sigZones_2$zone_id)

			# get overlaps of zones 
			ranges <- subsetByOverlaps(gr1, gr2)
			hits <- findOverlaps(gr1, gr2)
			rsid <- CharacterList(split(gr2$names[subjectHits(hits)], queryHits(hits)))
			mcols(ranges) <- DataFrame(mcols(ranges), rsid)
			ranges

			# convert zone ids to dt
			df <- data.table(data.frame(mcols(ranges)))
			setnames(df, c("names", "rsid"), c("trait_1_zone_id", "trait_2_zone_id"), skip_absent=TRUE)	

			# first remove list characters from string then split into long from wide
			df$trait_2_zone_id <- gsub('["c()"]', '', df$trait_2_zone_id )
			
			# if there are no matches, then just put NA for both, we can filter this out later
			if (nrow(df) > 0) {
			df <- df[, lapply(.SD, function(x) unlist(tstrsplit(x, ",", fixed=TRUE))), by = trait_1_zone_id][!is.na(trait_2_zone_id)]

			}else{df <- data.table(cbind(trait_1_zone_id = NA, trait_2_zone_id = NA))}

			# save overlap size also. query is sigzone 1 and subject is sigzone 2
			#df$overlap_size <- width(GenomicRanges::intersect(gr1[queryHits(hits)], gr2[subjectHits(hits)])) 

			#int_1 <- intersect(gr1[queryHits(hits)], gr2[subjectHits(hits)])
			#int_2 <- GenomicRanges::intersect(gr1, gr2)

			#df$overlap_size <- width(GenomicRanges::intersect(gr1, gr2))
			df$trait_1 <- trait
			df$trait_2 <- trait2

			# save the trait pairs and data frame
			name[counter] <- paste(trait, trait2, sep = "_")
			res_list[[counter]] <- df
			
			counter = counter + 1
	  	}


	}
	

	zone_overlap <- rbindlist(res_list)

	zone_overlap$treatment <- subset

	# this gives the long data frame
	#zone_overlap[, .N, by = c("trait_1", "trait_2")]

	# this gives the matrix
	#z <- dcast(zone_overlap, trait_1 ~ trait_2, value.var = "trait_1_zone_id")

	# combine both trait and id and columns and bind for full list
	r_1 <- data.table(cbind(id = zone_overlap$trait_1_zone_id, trait = zone_overlap$trait_1))
	r_2 <- data.table(cbind(id = zone_overlap$trait_2_zone_id, trait = zone_overlap$trait_2))

	r_full <- rbindlist(list(r_1, r_2))

	# remove duplicate by id and trait and na if needed
	total_df <- unique(r_full, by=c("id", "trait"))
	total_df <- na.omit(total_df)

	# get number of times id was shown, thus pleio
	total_df[, pleio := (.N), by = id]
	fwrite(total_df, sep = "\t", paste0(subset, "id_trait_pleio_table_sigzone_only.tsv"))

	# get pleio table
	y <- total_df[, (.N/pleio), by = pleio]
	fwrite(y, sep = "\t", paste0(subset, "_pleio_table_sigzone_only.tsv"))

	print(paste0("finished pleio files for subset ", subset, "- - - - - - - - - - - - - - - - - - - - - - - - -"))
}
