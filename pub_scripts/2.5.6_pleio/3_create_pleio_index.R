if (!require(tidyverse)) install.packages('tidyverse')
if (!require(vroom)) install.packages('vroom')
if (!require(data.table)) install.packages('data.table')
if (!require(ggplot2)) install.packages('ggplot2')
if (!require(viridis)) install.packages('viridis')
if (!require(ggpubr)) install.packages('ggpubr')
if (!require(igraph)) install.packages('igraph')

file_wd <- "~/mol_eco_2024/2.5.5_gwas/"
out_wd <- "~/mol_eco_2024/2.5.5_gwas/"

treat <- c("LHB", "LHH", "LNHB", "LNHH", "THB", "THH", "TNHB", "TNHH")

pleio_list <- vector("list", length = length(treat))

for(i in 1:length(treat)) {
	subset = treat[i]
	print(paste0("loading files for subset ", subset))

	# set wd to place outfile there
	out_dir <- paste(out_wd, subset, "/", sep = "")
	setwd(out_dir)
	#getwd()

	# read in datasets for plotting
	y <- fread(sep = "\t", paste0(subset, "_pleio_table_sigzone_only.tsv"))
	y[, treatment := subset]
	setkey(y, pleio)
	y = y[-1,]

	setnames(y, c("V1"), c("N"), skip_absent=TRUE)	

	y[, s := pleio*N]
	sum_pleio = sum(y$N)
	y[, weighted_s := s/sum_pleio]

	ind <- sum(y$weighted_s)#/length(y$pleio)

	pleio_list[[i]] <- data.table(ind)



}

pleio_index <- rbindlist(pleio_list)
pleio_index$treatment <- treat

setwd(out_wd)
fwrite(pleio_index, sep = ",", paste0("pleio_index_sigzone_only.csv")) 

t_mat_long <- data.table(melt(trait_matrix, value.name = "snps_shared", varnames=c('trait_1', 'trait_2')))