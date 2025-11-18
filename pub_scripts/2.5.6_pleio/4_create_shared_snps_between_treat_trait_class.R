# test if morpho traits share more SNPs with morpho traits
if (!require(tidyverse)) install.packages('tidyverse')
if (!require(vroom)) install.packages('vroom')
if (!require(data.table)) install.packages('data.table')
if (!require(igraph)) install.packages('igraph')
if (!require(bipartite)) install.packages('bipartite')


# loop through treatments
treat <- c("LHB", "LHH", "LNHB", "LNHH", "THB", "THH", "TNHB", "TNHH")

# set wd to place outfile there
file_wd <- "~/mol_eco_2024/2.5.5_gwas/"
out_wd <- "~/mol_eco_2024/2.5.5_gwas/"

# create list for matrices
mat_list <- vector("list", length = length(treat))

# pairwise loop to assess network, load in all treat networks into list and mds
for(i in 1:(length(treat))) {
	subset = treat[i]
	out_dir <- paste(out_wd, subset, "/", sep = "")
	setwd(out_dir)
	print(paste0("loading files for subset ", subset, "- - - - - - - - - - - - - - - - - - - - - - - - -"))

	# read in treatment file
	pleio_df <- fread(sep = "\t", paste0(subset, "id_trait_pleio_table_sigzone_only.tsv"))
	pleio_df <- unique(pleio_df, by = c('id', "trait"))

	### plot only where pleio tropic, ie 2 or more connections of snp to traits - - - - - - - - - - - - - - - - - - - - - - - - - 
	d_pleio <- pleio_df[pleio >= 2]
	pleio_weights <- unlist(as.vector(abs(d_pleio[, 3])))
	d_pleio <- d_pleio[,1:2]
	d_pleio_mat = table(d_pleio)
	class(d_pleio_mat) <- "matrix" # And we convert it from a table to a matrix

	# transpose to create the trait matrix
	trait_matrix = t(d_pleio_mat) %*% d_pleio_mat

	# melt to long, add treatment then add to list
	trait_matrix[lower.tri(trait_matrix, diag = TRUE)] <- NA
	t_mat_long <- data.table(melt(trait_matrix, value.name = "snps_shared", varnames=c('trait_1', 'trait_2')))
	t_mat_long$treat <- subset

	mat_list[[i]] <- t_mat_long
}

# create big boy shared snps data table
t_mat <- rbindlist(mat_list)
t_mat <- na.omit(t_mat)

# add trait class
trait_df <- fread(sep = "\t", paste0("~/mol_eco_2024/ch_1_sub/2.5.5_gwas/trait_class.tsv")) 

# consalidate all volatiles to volatile
trait_df$trait_class <- ifelse(trait_df$trait_class != "morphology", "volatile",trait_df$trait_class)
trait_df$trait_class <- as.factor(trait_df$trait_class)

# add trait class for trait 1
setkey(t_mat, "trait_1")
setkey(trait_df, "trait")
t_mat <- t_mat[trait_df]
setnames(t_mat, c("trait_class"), c("trait_class_1"))

# add trait class for trait 2
setkey(t_mat, "trait_2")
setkey(trait_df, "trait")
t_mat <- t_mat[trait_df]
setnames(t_mat, c("trait_class"), c("trait_class_2"))

# create the shared class col, change so that vol_morph is same as morph_vol, ie order doesnt matter
t_mat[, shared_class := paste(trait_class_1, trait_class_2, sep = "_")]
t_mat$shared_class <- ifelse(t_mat$shared_class == "volatile_morphology", "morphology_volatile", t_mat$shared_class)

t_mat <- na.omit(t_mat)


setwd(out_wd)
fwrite(t_mat, sep = "\t", paste0("shared_snps_between_traits_sigzone_only.tsv"))