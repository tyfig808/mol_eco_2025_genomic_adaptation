### create emmax files for gwas
setwd("~/mol_eco_2024/2.5.5_gwas")

if (!require(tidyverse)) install.packages('tidyverse')
if (!require(data.table)) install.packages('data.table')


### load in phenotypic data, run pca on traits proven to be sig
d <- fread("emmax_phenotypic_data_thomas.csv")

### here are the traits needed for the pheno file
morphology <- c("Heightd22","Heightd34","Leaf-size","Floweringtime","Brancheslength","Nbbranches","Flowerperdaysincefloweringtime","NbFlowers",
        "Flowerdiameter","Petal-length","Petal-width","Sepale-length","Style","Stamen") # missing nectar amount

mating_sys <- c() # couldnt find traits in dataset
#b_assays <- c("Visits") # because g1 didnt have this, causing to remove g1 from data, best to remove
volatiles <- c("LNBenzaldehyde","LNMethyl-benzoate","LNPhenylethyl-Alcohol","LNPhenylacetaldehyde","LNp-Anisaldehyde","LNMethyl-salicylate")
terp <- c("LNa-Farnesene-E")
fatty <- c("LN3-Hexen-1-ol-acetate-Z")
nitro <- c("LNMethyl-anthranilate","LNBenzyl-nitrile","LNIndole","LN2-Amino-benzaldehyde")
sulfer <- c("LN1-butene-4-isothiocyanate")

# split into the cols we need
d<- as.data.frame(d)
x <- d[, c(morphology, volatiles, terp, fatty, nitro, sulfer)]
x <- cbind(x, ind_id = d$PlantlabelDNA)

# move id up to first col
x<- x[,c(ncol(x),1:(ncol(x)-1))]

# change back to write
x <- data.table(x)

# change pheno of 27A to 27A_c the resequneced sample
x<- x %>% 
        mutate(ind_id = recode(ind_id, RA_TL_H_B_27A = "RA_TL_H_B_27")) %>%
        mutate(ind_id = recode(ind_id, RA_TL_H_B_27A_c = "RA_TL_H_B_27A")) %>%
        mutate(ind_id = recode(ind_id, RA_TL_H_B_27 = "RA_TL_H_B_27A_c")) 

fwrite(x, sep = "\t", paste0("emmax_pheno_data.tsv"))

# create covariate files and check length
phen <- x

allOutputFiles <- list.files(pattern = "G1.numchr.tfam", recursive = TRUE, full.names = TRUE)
f_names <- str_sub(allOutputFiles,9,-27)
str_sub(f_names,5,5) <- "" 

k_names <- str_sub(allOutputFiles,3,-6)

for(i in (1:length(allOutputFiles))) {
                # redo name of this sample as it will give problems, changed back after merge
                print(paste0("subset file: ", f_names[i], " being processed"))
                
                #setnames(phen, c("LN3_Hexen_1_ol_acetate_(Z)_"), c("LN3_Hexen_1_ol_acetate_Z"), skip_absent=TRUE)

                # read in samples file with which we want to merge pheno data, remove the a or b at the end and remove single obs
                dt <- fread(allOutputFiles[i], header = FALSE) 
                dt[,id := str_sub(dt$V2,1,-2)]
                
                # for thb, becasue there is the resequenced 
                #dt <- fread("/shares/schiestl.systbot.uzh/variants/data/tfig/gzp_vcf/gwas_files/subsetTHBAsubsetTHBB_G1_original.numchr.tfam", header = FALSE) 
                #dt[,id := str_sub(dt$V2,1,-2)]

                # filter and check if data tables are the same length
                sub <- phen[dt, on = .(ind_id == id)]

                y = f_names[i]
                ifelse(nrow(dt) == nrow(sub), print(paste0(y, ": equal number of rows: ", nrow(dt)," = ", nrow(sub))), 
                                print(paste0(y, " !!! wrong number of rows !!! sample file = ", nrow(dt)," | phen file = ", nrow(sub))))
                
                # if need to add family id
                #c <- cbind(c, gen = str_split_i(c$ind_id, "_", 2))


                # make it so the sub pheno data table has the names of the genetic samples file, 
                sub[, c("ind_id", "V1", "V3", "V4", "V5", "V6"):=NULL]
                setnames(sub, c("V2"), c("ind_id"), skip_absent=TRUE)
                setcolorder(sub, neworder = "ind_id")

                # for THB, there is an original file, load instead, return the pheno file back to the loop and write pheno file for gwas
                #sub<- sub %>% 
                #       mutate(ind_id = recode(ind_id, RA_TL_H_B_27A = "RA_TL_H_B_27")) %>%
                #       mutate(ind_id = recode(ind_id, RA_TL_H_B_27A_c = "RA_TL_H_B_27A")) %>%
                #       mutate(ind_id = recode(ind_id, RA_TL_H_B_27 = "RA_TL_H_B_27A_c")) 
                #sub <- sub[ind_id != "RA_TL_H_B_27A"]


                fwrite(sub, sep = "\t", paste0(y,"_phen_emmax.txt"))

                # create kinship file with emmax
                system(paste0("~/mol_eco_2024/emmax/emmax-kin-intel64 -v -s -d 10 ~/mol_eco_2024/", k_names[i]))
}


