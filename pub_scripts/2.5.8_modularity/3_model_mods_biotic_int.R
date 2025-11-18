# run simulations of distance matrix to see if we can say two matrixses are dissimilar
if (!require(tidyverse)) install.packages('tidyverse')
if (!require(vroom)) install.packages('vroom')
if (!require(data.table)) install.packages('data.table')
if (!require(car)) install.packages('car')
if (!require(glmmTMB)) install.packages('glmmTMB')
if (!require(DHARMa)) install.packages('DHARMa')
if (!require(fitdistrplus)) install.packages('fitdistrplus')

# loop through treatments
treat <- c("LHB", "LHH", "LNHB", "LNHH", "THB", "THH", "TNHB", "TNHH")

# set wd to place outfile there
file_wd <- "~/mol_eco_2024/2.5.5_gwas/"
out_wd <- "~/mol_eco_2024/2.5.5_gwas/"

# create list to add the data table into, we will rbind after
m_list <- vector("list", length = length(treat))

# retrive modules from 3 methods, loop over treatments, has problems with ubuntu for whatever reason, glpk or something 
for(i in 1:(length(treat))) {

	subset = treat[i]
	out_dir <- paste(out_wd, subset, "/", sep = "")
	setwd(out_dir)
	print(paste0("loading files for subset ", subset, "- - - - - - - - - - - - - - - - - - - - - - - - -"))

	# load in shared modules, get treatment and length
	s_mod <- fread(sep = "\t", paste0(subset, "_shared_modules_sigzone_only.tsv"))

	m_list[[i]] <- data.table(nrow(s_mod))
}

# bind to list add treat
n_mod <- rbindlist(m_list)
n_mod$treat <- treat

# add factor cols
n_mod$soil <- str_sub(n_mod$treat, 1, 1)
n_mod$herb <- str_sub(n_mod$treat, 2, -2)
n_mod$bee <- str_sub(n_mod$treat, -1)

# add biotic cols
n_mod$h_bio <- ifelse(n_mod$herb == "H", 1, 0) 
n_mod$b_bio <- ifelse(n_mod$bee == "B", 1, 0)
n_mod$s_bio <- rowSums(n_mod[,(length(n_mod)-1):(length(n_mod))])

# check whether biotic affect number of mods with soil int
descdist(n_mod$V1, discrete = TRUE)
m_mod <- glmmTMB(V1 ~ s_bio * soil , data = n_mod)
Anova(m_mod)
#Response: V1
#            Chisq Df Pr(>Chisq)  
#s_bio      1.7561  1     0.1851  
#soil       0.8780  1     0.3487  
#s_bio:soil 4.8780  1     0.0272 *

simulateResiduals(fittedModel = m_mod, plot = T, n = 100)


# check relationship with n mods and size
setwd(out_wd)
size_dt <- fread(sep = "\t", paste0("network_size_treatments_sigzone_only.tsv"))



size_mod <- glmmTMB(1/size_dt$size ~ n_mod$V1, family = beta_family())
Anova(size_mod)
#          Chisq Df Pr(>Chisq)  
#n_mod$V1 2.8186  1    0.09318 .

simulateResiduals(fittedModel = size_mod, plot = T, n = 100)

