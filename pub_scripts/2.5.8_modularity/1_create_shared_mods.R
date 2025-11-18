# run simulations of distance matrix to see if we can say two matrixses are dissimilar
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

# create color scale for trait class
trait_df <- fread(sep = "\t", paste0("~/mol_eco_2024/2.5.5_gwas/trait_class.tsv")) 
trait_df$trait_class <- as.factor(trait_df$trait_class)

col_max <- length(unique(trait_df$trait_class))
col <- data.table(color = viridis_pal(alpha = 1, option = "turbo")(col_max))
col = cbind(col, trait_class = unique(trait_df$trait_class))
setkey(col, "trait_class")
setkey(trait_df, "trait_class")
trait_df <- col[trait_df]

# load big boy function
source('~/function_shared_mods_sigzone_only.R')

# retrive modules from 3 methods, loop over treatments, has problems with ubuntu for whatever reason, glpk or something 
for(i in 1:(length(treat))) {

subset = treat[i]
# run function to extract the shared modules
final_res <- shared_mods(i)
#final_res <- final_res[order(-len)]

# give new subset name so we can collect them all later
assign(paste0("final_res_", subset), final_res)

# write the file before the next part of the loop 
fwrite(final_res, sep = "\t", paste0(subset, "_shared_modules_sigzone_only.tsv"))

}