# create full fst sep within reps, ie only G1_A against LHB_A
if (!require(tidyverse)) install.packages('tidyverse')
if (!require(vroom)) install.packages('vroom')
if (!require(data.table)) install.packages('data.table')


# set wd read in files
setwd("/shares/schiestl.systbot.uzh/variants/data/tfig/ch_1_sub/2.5.1_fst_pair")

# get files and create names, 
fst_files <- list.files(pattern = "_pop.windowed.weir.fst", recursive = TRUE, full.names = TRUE)

# get rep vector
rep_1 <- str_split_i(fst_files, "_", 1)
rep_1 <- str_split_i(rep_1, "/", 2)

rep_2 <- str_split_i(fst_files, "_", 3)

# get pop vectors
pop_1 <- str_split_i(fst_files, "_", 2)
pop_2 <- str_split_i(fst_files, "_", 4)

# remove the second soil from but G1
g_1 <- which(pop_1!="G1")
str_sub(pop_1[g_1],2,2) <- ""  # remove the second L from all groups but G1

g_2 <- which(pop_2!="G1")
str_sub(pop_2[g_2],2,2) <- "" # remove the second L from all groups (there is no G1)

# loop over rep a and rep b to create seperature full file
r <- unique(c(rep_1, rep_2))

for(l in 1:length(r)) {
       # load rep name, create subset where both pairwise have the same rep, and either rep A or B, whatver the loop says
       r_sub <- r[l]
       print(paste0("loading files for rep: ", r_sub))
       
       a_a <- which(rep_1 == rep_2 & rep_1 == r_sub)

       # create a list to save the rep files and write out later
       fst_list <- vector("list", length = length(a_a))

       # loop over the subset of rep files
       for(i in 1:length(a_a)) {
              # run over the index
              j = a_a[i]
              print(paste0("loading files for pop 1: ", pop_1[j], " and for pop 2: ", pop_2[j]))

              # read in fst, add colums for treatments
              fst <- fread(sep = "\t", fst_files[j])

              # add pop 1 and 2 cols with the factors to run the anova on 
              fst$pop_1 <- pop_1[j]
              fst$pop_2 <- pop_2[j] 
              
              # add the other factors, need to think a bit about how to structure
              fst$gen_1 <- ifelse(pop_1[j]=="G1", 1, 10)
              fst$soil_1 <- ifelse(pop_1[j]=="G1", "NA", str_sub(pop_1[j],1,1))
              fst$herb_1 <- ifelse(pop_1[j]=="G1", "NA", str_sub(pop_1[j],2,-2))
              fst$bee_1 <- ifelse(pop_1[j]=="G1", "NA", str_sub(pop_1[j],-1,-1))

              # repeat for pop 2
              fst$gen_2 <- ifelse(pop_2[j]=="G1", 1, 10)
              fst$soil_2 <- ifelse(pop_2[j]=="G1", "NA", str_sub(pop_2[j],1,1))
              fst$herb_2 <- ifelse(pop_2[j]=="G1", "NA", str_sub(pop_2[j],2,-2))
              fst$bee_2 <- ifelse(pop_2[j]=="G1", "NA", str_sub(pop_2[j],-1,-1))

              # create col with both combined
              fst$gen_x <- paste(fst$gen_1, fst$gen_2 , sep = "_")
              fst$soil_x <- ifelse(pop_1[j]=="G1", paste("G1", fst$soil_2 , sep = "_"), paste(fst$soil_1, fst$soil_2 , sep = "_"))
              fst$herb_x <- ifelse(pop_1[j]=="G1", paste("G1", fst$herb_2 , sep = "_"), paste(fst$herb_1, fst$herb_2 , sep = "_"))
              fst$bee_x <- ifelse(pop_1[j]=="G1", paste("G1", fst$bee_2 , sep = "_"), paste(fst$bee_1, fst$bee_2 , sep = "_"))

              fst_list[[j]] <- fst

       }


       full_fst <- rbindlist(fst_list) 

       # remove scaffolds from fst dataset
       chr <- c("A01", "A02", "A03", "A04", "A05",
              "A06", "A07", "A08", "A09", "A10") 

       full_fst <- full_fst[full_fst$CHROM %in% chr, ] # take all that match the chr object

       # write out
       fwrite(full_fst, sep = "\t", paste0(r_sub, "_fst_full_data.tsv"))

}


rep <- c("RA", "RB")

# create list to rbind
fst_list <- vector("list", length = length(rep))

# loop over reps, read in and the bind 
for(l in 1:length(rep)) {
       r_sub <- rep[l]
       print(paste0("loading files for rep: ", r_sub))
       
       # read in fst, add rep column and add to list
       rep_fst <- fread(sep = "\t", paste0(r_sub, "_fst_full_data.tsv"))
       rep_fst$rep <- r_sub

       fst_list[[l]] <- rep_fst

}

full_fst <- rbindlist(fst_list)

fwrite(full_fst, sep = "\t", paste0("total_rep_sep_fst_full_data.tsv"))
