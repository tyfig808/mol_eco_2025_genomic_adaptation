# create sup mod table, mark with x if in there and 0 if not
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

# load in trait df to get traits
trait_df <- fread(sep = "\t", paste0("~/mol_eco_2024/2.5.5_gwas/trait_class.tsv")) 
t <- trait_df$trait

# create list to add the data table into, we will rbind after
t_list <- vector("list", length = length(treat))
len_list <- vector("list", length = length(treat))

# retrive modules from 3 methods, loop over treatments, has problems with ubuntu for whatever reason, glpk or something 
for(i in 1:(length(treat))) {

subset = treat[i]
out_dir <- paste(out_wd, subset, "/", sep = "")
setwd(out_dir)
print(paste0("loading files for subset ", subset, "- - - - - - - - - - - - - - - - - - - - - - - - -"))

# load in shared modules, get treatment and length
res <- fread(sep = "\t", paste0(subset, "_shared_modules_sigzone_only.tsv"))

t_list[[i]] <- res

}


t <- rbindlist(t_list)

setwd(out_wd)
fwrite(t, sep = ",", paste0("sup_mod_list_sigzone_only.csv")) 
