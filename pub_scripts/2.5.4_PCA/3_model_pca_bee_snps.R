### see if number of snps relates to PC1 and PC2 combined, maybe get the vector from G1?
if (!require(tidyverse)) install.packages('tidyverse')
if (!require(vroom)) install.packages('vroom')
if (!require(data.table)) install.packages('data.table')
if (!require(glmmTMB)) install.packages('glmmTMB')
if (!require(DHARMa)) install.packages('DHARMa')
if (!require(car)) install.packages('car')

# read in cmh snp pca
setwd("~/mol_eco_2024/2.5.4_pca")

#dt <- fread(paste0("pca_final_fdr_snp_no_indel.eigenvec"))

full_pca <- fread(paste0("pca_final_total_snp_no_indel.eigenvec")) 
sig_pca <- fread(paste0("pca_clear_cmh_sig_no_indel.eigenvec"))


pca_list <- list(full_pca, sig_pca)
pca_names <- c("full_pca", "sig_pca")

cen_list <- vector("list", length = length(pca_list))

# for the number of PCAs calcuate the distance from G1 from each treatment, do it within Replicates
for(i in 1:(length(pca_list))) {
  dt <- pca_list[[i]]
  setnames(dt, "#IID", "id", skip_absent=TRUE)

  ### replace 27A with resequence 
  dt <- dt %>%
    #subset(sigzone_pca$V1 != "RA_TL_H_B_27A") %>%
    mutate(id = dplyr::recode(id, "RA_TL_H_B_27A_c" = "RA_TL_H_B_27A"))

    # get rep, then soil, remove the second soil if it is not G1
  dt$rep<- str_split_i(dt$id, "_", 1) 
  dt$soil<- str_split_i(dt$id, "_", 2)

  g_1 <- which(dt$soil!="G1")
  str_sub(dt$soil[g_1],2,2) <- ""  # remove the second L from all groups but G1

  dt$herb<- ifelse(dt$soil != "G1", str_split_i(dt$id, "_", 3), NA)
  dt$bee<- ifelse(dt$soil != "G1", str_split_i(dt$id, "_", 4), NA)

  # add treatment and group for plotting later
  dt[, treatment := paste(soil, herb, bee, sep = '_')]
  dt$treatment<- ifelse(dt$treatment == "G1_NA_NA", "G1", dt$treatment)

  dt[, group := paste(rep, treatment, sep = '_')]

  cen_sig <- data.table(aggregate(cbind(dt$PC1,dt$PC2)~group,dt,mean))
  # add rep col, then perfom within rep
  cen_sig$rep<- str_split_i(cen_sig$group, "_", 1) 

  # add generation to do operation on 
  cen_sig$gen <- str_split_i(cen_sig$group, "_", 2)  
  cen_sig$gen <- ifelse(cen_sig$gen != "G1", "G10", cen_sig$gen) 

  # loop over rep to get measure
  rep <- unique(cen_sig$rep)
  sig_list <- vector("list", length = length(rep))

  for(j in 1:(length(rep))) {
    r <- rep[j]
    
    # for this rep, use this G1 and G10
    b <- cen_sig[rep == r & gen == "G1"]
    b_10 <-cen_sig[rep == r & gen == "G10"]

    # Take the euclidean distance, ie the hypotenuse
    b_10[, evoL_dist := sqrt(((V1-b$V1)**2) +((V2-b$V2)**2))]

    # add to list so can combine with rep
    sig_list[[j]] <- b_10

    }

  # take average across rep for treatment
  b_10 <- rbindlist(sig_list)
  b_10$treat <- str_sub(b_10$group,4,-1)
  treat_dt <- b_10[, mean(evoL_dist), by = treat]

  # add treatment columns and which dataset it came from
  treat_dt$soil <- str_split_i(treat_dt$treat, "_", 1)
  treat_dt$herb <- str_split_i(treat_dt$treat, "_", 2)
  treat_dt$bee <- str_split_i(treat_dt$treat, "_", 3)

  treat_dt$data_set <- str_split_i(pca_names[i], "_", 1) 

  cen_list[[i]] <- treat_dt
}

# bind to full list
evo_change <- rbindlist(cen_list)

# check only for sig
treat_dt <- evo_change[data_set == "sig"]

# run model, bees treatments diverged more from G1 
bee_snp_mod <- glmmTMB(V1 ~ bee, treat_dt)
Anova(bee_snp_mod)
#     Chisq Df Pr(>Chisq)    
#bee 10.456  1   0.001223 **
