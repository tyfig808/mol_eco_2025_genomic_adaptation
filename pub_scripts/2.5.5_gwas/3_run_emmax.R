###  array for emmax over treatments
setwd("~/mol_eco_2024/2.5.5_gwas")

if (!require(tidyverse)) install.packages('tidyverse')
if (!require(vroom)) install.packages('vroom')
if (!require(data.table)) install.packages('data.table')

### for emmax we need phen and kin files, and then a tped tfam
phen_files <- list.files(pattern = "_phen_emmax.tsv", recursive = TRUE, full.names = TRUE)
kin_files <- list.files(pattern = "G1.numchr.aIBS.kinf", recursive = TRUE, full.names = TRUE)

## for tfam files, need to create emmax pheno file, although the emmax wants the tped prefix, same number of characters so it works
t_files <- list.files(pattern = "_G1.numchr.tfam", recursive = TRUE, full.names = TRUE)

f_names <- str_split_i(phen_files, "_", 1)
f_names <- str_sub(f_names,3,6)

# for tfam files, emmax just wants the prefix
t_names <- str_split_i(t_files, "/", 2)
t_names<- str_sub(t_names,1,-6)

k_names <- str_split_i(kin_files, "/", 2)

# use this to run with array in bash at same time
args <- commandArgs(trailingOnly=T)
j <- eval(parse(text=args[1]))
#j = 5

print(paste0("loading files for subset ", f_names[j]))

# load in phen
phen <-  fread(phen_files[j]) # think we need the head here
#phen=read.table(phen_files[j], h=T,na.strings=".")
phen <- as.data.frame(phen)

# set wd to place outfile there
file_wd <- "~/mol_eco_2024/2.5.5_gwas"
out_wd <- "~/mol_eco_2024/2.5.5_gwas"
out_dir <- paste(out_wd, f_names[j], "/", sep = "")
setwd(out_dir)


### emmax over trait
#trait <- "Heightd30"
for(trait in colnames(phen)[2:ncol(phen)]) {

	# load tfam file and add the phenotypic data each time, then write it, emmax then loads it in
	X <-  fread(paste("/", t_files[j], sep = ""), select = c(1,2))

	# subset phen, merge and write, can reduce to one line with just cbind, but this should be safer as it matches the names
	phen_t <- phen[, c("ind_id", trait)]
	X <- X[phen_t[,c("ind_id",trait)], on = .(V2 == ind_id)]
	X <- as.data.frame(X)
	write.table(X, paste(out_dir, f_names[j], "_", trait, ".phen.txt", sep=""), quote=F, row.names=F, col.names=F, sep=" ")
	#fwrite(X, paste(out_dir, f_names[j], "_", trait, ".phen", sep=""), row.names=FALSE, col.names=FALSE, sep=" ")

	print(paste0("Running gwas for subset ", f_names[j], " and for trait: ", trait))
	system(paste("~/mol_eco_2024/emmax/emmax-intel64 -v -d 10",
		" -t ", file_wd, t_names[j], 
		" -p ", out_dir, f_names[j], "_", trait, ".phen.txt", 
		" -k ", file_wd, kin_files[j], 
		" -o ", out_dir, f_names[j], "_", trait, ".emmax", sep=""))

} 


#% emmax -v -d 10 -t [tped_prefix] -p [pheno_file] -k [kin_file] -o [out_prefix]




