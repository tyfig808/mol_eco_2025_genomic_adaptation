# subset full freq file into cmh subset, for fdr and updated ne
if (!require(tidyverse)) install.packages('tidyverse')
if (!require(data.table)) install.packages('data.table')

# set some wd
file_wd <- "/shares/schiestl.systbot.uzh/variants/data/tfig/gzp_vcf/gwas_files/"
out_wd <- "/shares/schiestl.systbot.uzh/variants/data/tfig/gwas/"

# this is the kb range used
r <- 1000
kb <- r/1000

# loop through treatments
treat <- c("LHB", "LHH", "LNHB", "LNHH", "THB", "THH", "TNHB", "TNHH")
pruned <- paste(treat, "A", sep = "")

# create table to convert cmh chrome to just 1 or 2
x1 <- c(1,2,3,4,5,6,7,8,9,10)
x2 <- c("A01","A02","A03","A04","A05","A06","A07","A08","A09","A10")
look <- data.table(x1,x2)
setkey(look, x2)

# cols to take for af
myVector <- c("chr", "pos", "i.fdr", "G10a_af", "G1a_af", "id")

# create list to save to
treat_list <- vector("list", length = length(treat))

# loop over treatment and load in all cmh
for(i in 1:(length(treat))) {
	subset = treat[i]
	print(paste0("loading files for subset ", subset, "- - - - - - - - - - - - - - - - - - - - - - - - -"))

	# read in snps sig in cmh with updated ne
	setwd("/shares/schiestl.systbot.uzh/variants/data/tfig/freq_files")
	cmh <- fread(sep = "\t", paste0("subset_", subset, "cmh_fdr_sig_updated_ne.tsv"))
	setkey(cmh, ID)

	# read in list of from plink that are present after ld analyses
	setwd("/shares/schiestl.systbot.uzh/variants/data/tfig/gzp_vcf")
	p <- pruned[i]
	clean_snps <- fread(sep = "\t", header = FALSE, paste0("pruned_", p, ".prune.in"))

	# for in file sep by :, here we change to _
	clean_snps$V1 <- sub(':','_', clean_snps$V1)
	setkey(clean_snps, V1)

	# merge and take where match in both, add to cmh_list for merging later
	cleaned_cmh <- cmh[clean_snps, nomatch = NULL]
	setkey(cleaned_cmh, chr)

	# get numerical chrome for sigzone matching
	cleaned_cmh <- cleaned_cmh[look]
	setnames(cleaned_cmh, c("x1"), c("chr_num"), skip_absent=TRUE)	
	setkey(cleaned_cmh, chr_num, pos)

	# moving into gwas subset directories
    out_dir <- paste(out_wd, subset, "/", sep = "")
    setwd(out_dir)

    ### get sigzone files with beg and end and full local score with distances
    kb_files <- list.files(pattern = "kb_range.tsv", recursive = TRUE, full.names = TRUE)
    trait_names <- str_split_i(kb_files, "_", 2)

    # keep track of traits and for treatment totals
    trait_list <- vector("list", length = length(trait_names))
    
    # trait loop
    for(k in 1:length(kb_files)) {

        trait = trait_names[k]
        print(paste0("loading files for subset ", subset, " and for trait: ", trait))

        # load in snps in sigzones
        x <- fread(sep = "\t", paste0(subset, "_", trait, "_emmax_sigzones_", kb, "kb_range.tsv"))
        setkey(x, chr, pos)

        # get shared SNPs from both cmh and gwas
        shared <- x[cleaned_cmh, nomatch = NULL]
        shared <- shared[, ..myVector]
        setnames(shared, c("i.fdr"), c("cmh_fdr"), skip_absent=TRUE)	

        # if there are no snps in cmh and gwas then create empty row
        if (nrow(shared) == 0) {
        	shared <- data.table(t(data.table(c(NA,NA,NA,NA,NA,NA))))
        	setnames(shared, c("V1", "V2", "V3", "V4", "V5", "V6"), 
        		c("chr", "pos", "cmh_fdr", "G10a_af", "G1a_af", "G10b_af", "G1b_af", "id"), skip_absent=TRUE)
        }
        
        # add the trait and add to list 
        shared$trait <- trait
        trait_list[[k]] <- shared
    }

    # after running all traits for the treatment, bind to list and add treatment col, add to treatment list
    shared_trait <- rbindlist(trait_list)
    shared_trait$treatment <- subset
    treat_list[[i]] <- shared_trait 

}

# after all treatments, bind to list, get full 
full_shared <- rbindlist(treat_list)

# add af columns for transformed difference
full_shared[, a_diff := abs(G10a_af-G1a_af)]
full_shared[, b_diff := abs(G10b_af-G1b_af)]
full_shared[, mean_diff := mean(a_diff, b_diff)]
full_shared[, trans_af_diff := (2.0*asin(sqrt(mean_diff)))]

setwd(out_wd)
fwrite(full_shared, sep = "\t", paste0("shared_cmh_af_all_treatments_traits.tsv")) 





    