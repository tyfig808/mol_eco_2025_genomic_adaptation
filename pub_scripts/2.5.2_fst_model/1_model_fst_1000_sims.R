# mixed model on fst to say which factor matters the most, boxcox works well, rrpp model just run once
if (!require(tidyverse)) install.packages('tidyverse')
if (!require(data.table)) install.packages('data.table')
#if (!require(DHARMa)) install.packages('DHARMa')
if (!require(car)) install.packages('car')
if (!require(RRPP)) install.packages('RRPP')
if (!require(MASS)) install.packages('MASS')


# set wd
setwd("~/mol_eco_2024/2.5.1_fst_pair")

# read in full fst, this one has rep sep
full_fst <- fread(sep = "\t", paste0("total_rep_sep_fst_full_data.tsv"), stringsAsFactors=TRUE)

# take only G10
g10_fst <- full_fst[pop_1 != "G1"]

rm(full_fst)

# do boxcox for G1 and use 
g10_filt <- g10_fst[WEIGHTED_FST>0]
bc <- boxcox(lm(g10_filt$WEIGHTED_FST ~ 1))
lambda <- bc$x[which.max(bc$y)]

rm(g10_filt)
gc()

# create pair treatment file to get number of obs per treat
g10_fst$pair_treat <- paste(g10_fst$pop_1, g10_fst$pop_2, sep = "_")

g10_fst$soil_x <- droplevels(g10_fst$soil_x)
g10_fst$herb_x <- droplevels(g10_fst$herb_x)
g10_fst$bee_x <- droplevels(g10_fst$bee_x)
g10_fst[, trans_fst :=((WEIGHTED_FST^lambda-1)/lambda)]

# get number of obs per random sample
total_obs <- nrow(g10_fst)/length(unique(g10_fst$pair_treat))/length(unique(g10_fst$rep)) # number of rows divided by treatments pairs divided by rep size
per <- .50
n <- 1000

# set sim number and create data table to save 
nsim <- 1000
f_list_trans <- vector("list", length = nsim)
f_list <- vector("list", length = nsim)
coef_list <- vector("list", length = nsim)
set.seed(123)

# save these coefficients
coef_vec <- c("soilL_T", "soilT_T", "herbH_NH", "herbNH_H", "herbNH_NH", "beeB_H", "beeH_B", "beeH_H",
	"soilL_T:herbH_NH", "soilT_T:herbH_NH", "soilL_T:herbNH_NH", "soilT_T:herbNH_NH", 
	"soilL_T:beeB_H", "soilT_T:beeB_H", "soilL_T:beeH_B", "soilT_T:beeH_B", "soilL_T:beeH_H", "soilT_T:beeH_H",
	"herbH_NH:beeB_H", "herbNH_H:beeB_H", "herbNH_NH:beeB_H", "herbH_NH:beeH_B", "herbNH_H:beeH_B", "herbNH_NH:beeH_B", 
	"herbH_NH:beeH_H", "herbNH_H:beeH_H",   "herbNH_NH:beeH_H")

for (i in 1:nsim) {
	# randomize sample, run model a few times, get dist of f values of interactions, we want two ways, and three ways
	sample_dt <- g10_fst[g10_fst[, .I[sample(.N, n)], by = pair_treat][[2]]] 
	sample_dt <- sample_dt[WEIGHTED_FST>0]

	# create rrpp data frame
	fst.rrpp <- rrpp.data.frame(
			fst = sample_dt$WEIGHTED_FST,
			trans_fst = sample_dt$trans_fst,
			soil = sample_dt$soil_x,
			herb = sample_dt$herb_x,
			bee = sample_dt$bee_x,
			rep = sample_dt$rep,
			chr = sample_dt$CHROM,
			window = sample_dt$BIN_START
	)

	fst_fit_trans <- lm.rrpp(trans_fst ~ soil * herb * bee 
			+ rep + chr/window, # random effects
				SS.type = "III", data = fst.rrpp, print.progress = FALSE, iter = nsim)

	# this would updates model to mixed and makes these terms random effects, we leave out because it just makes the F 1,
	# lets do this after because doesnt work in bash for some reason
	#trans_a <- anova.lm.rrpp(fst_fit_trans, effect.type = "F", error = er)

	# make anova table and save F values
	trans_a <- anova.lm.rrpp(fst_fit_trans, effect.type = "F") 
	trans_a_dt <- data.table(trans_a$table)
	f_list_trans[[i]] <- transpose(data.table(trans_a_dt$F))

	# save coeficcients
	c <- data.table(t(coef(fst_fit_trans)))
	coef_list[[i]] <- c[, ..coef_vec]


	if(i%%100==0){print(i)}


}

# create col vectors to drop
col_vec <- c("Residuals", "Total")

setwd("~/mol_eco_2024/2.5.2_fst_model")

# bind back to list for transformed with boxcox
f_val_trans <- rbindlist(f_list_trans)
colnames(f_val_trans) <- rownames(trans_a$table)
f_val_trans[, (col_vec) := NULL]
fwrite(f_val_trans, sep = "\t", paste0("fst_rrpp_perm_boxcox_g10_only_trans.tsv"))

# write out coefficients
coef_dt <- rbindlist(coef_list)
fwrite(coef_dt, sep = "\t", paste0("fst_rrpp_perm_boxcox_g10_only_trans_coef.tsv"))

# check pairs
#PWT <- pairwise(fst_fit, groups = interaction(fst.rrpp$soil, fst.rrpp$bee))
#summary(PWT, confidence = 0.95)