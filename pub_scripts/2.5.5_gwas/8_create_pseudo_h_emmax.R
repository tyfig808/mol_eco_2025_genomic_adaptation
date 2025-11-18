# get gene annotations for gwas traits
if (!require(data.table)) install.packages('data.table')
if (!require(stringr)) install.packages('stringr')

# loop over treatments 
treat <- c("LHB", "LHH", "LNHB", "LNHH", "THB", "THH", "TNHB", "TNHH")

# create treatment wd
out_wd <- "~/mol_eco_2024/2.5.5_gwas/"

# create treatment list to be saved into
treat_list <- vector("list", length = length(treat))

# loop over treatments
for(i in 1:length(treat)) {
    subset = treat[i]
    print(paste0("loading files for subset ", subset))

    # moving into gwas subset directories
    out_dir <- paste(out_wd, subset, "/", sep = "")
    setwd(out_dir)

    ### get sigzone files 
    sz_emmax <- list.files(pattern = ".emmax.reml", recursive = TRUE, full.names = TRUE)  
    trait_names <- str_split_i(sz_emmax, "_", 2)
    trait_names <- str_split_i(trait_names, '.emmax', 1)

    # create trait list and loop over traits
    trait_list <- vector("list", length = length(trait_names))

    # trait loop
    for(k in 1:length(trait_names)) {
        #k = 1

        trait = trait_names[k]
        print(paste0("loading files for subset ", subset, " and for trait: ", trait))

        # load in sigzones and create id col
        dt <- fread(sz_emmax[k])

        r <- data.table(t(dt))

        #r$treat <- subset
        r$trait <- trait

        trait_list[[k]] <- r

    }

    res <- rbindlist(trait_list)
    res$treat <- treat[i]

    treat_list[[i]] <- res

}

total_dt <- rbindlist(treat_list)

setnames(total_dt, c("V1", "V2", "V3", "V4", "V5", "V6"), 
    c("w_var", "wo_var", "va_ve", "va", "ve", "pseud_h"), skip_absent=TRUE)   

# create trait class col
setwd(out_wd)
trait_df <- fread(sep = "\t", paste0("trait_class.tsv")) 

# consalidate all volatiles to volatile
trait_df$trait_class <- ifelse(trait_df$trait_class != "morphology", "volatile",trait_df$trait_class)
trait_df$trait_class <- as.factor(trait_df$trait_class)

# add trait class for trait 1
setkey(total_dt, "trait")
setkey(trait_df, "trait")
total_dt <- total_dt[trait_df]


fwrite(total_dt, sep = "\t", paste0("emmax_estimates_all_treat_trait.tsv"))