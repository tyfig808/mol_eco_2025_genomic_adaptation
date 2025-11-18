# create haldane per trait file
setwd("~/mol_eco_2024/2.4_pheno_evo")

if (!require(tidyverse)) install.packages('tidyverse')
if (!require(data.table)) install.packages('data.table')
if (!require(nlme)) install.packages('nlme')

# read in phen
x <- fread(sep = "\t", paste0("emmax_pheno_data.tsv"))

# so these need to be logged, the volatile already logged
morphology <- c("Heightd22","Heightd34","Leaf-size","Floweringtime","Brancheslength","Nbbranches","Flowerperdaysincefloweringtime","NbFlowers",
        "Flowerdiameter","Petal-length","Petal-width","Sepale-length","Style","Stamen") # missing nectar amount

# natural log to remove hetero scad, add 1 to make log(0) = 0
x[, (morphology) := lapply(.SD, function(d) 
         (log(d+1)) ), .SDcols = morphology]


# first add columns to seperate
x$rep <- str_split_i(x$ind_id, "_", 1)

# first soil is evolved line = soil
x$soil<- str_split_i(x$ind_id, "_", 2)

# second soil is the one grown in = soil_plastic
x$soil_plastic<- ifelse(x$soil == "G1", str_split_i(x$ind_id, "_", 3), str_sub(x$soil,2,2))

# create soil_plastic, then remove from soil column 
str_sub(x$soil,2,2) <- ""
x$soil<- ifelse(x$soil == "G", "G1", x$soil)

# get herb and bee, create treatment again
x$herb<- str_split_i(x$ind_id, "_", 3)
x$herb <- ifelse(x$soil == "G1", NA, x$herb)

x$bee<- str_split_i(x$ind_id, "_", 4)
x$bee <- ifelse(x$soil == "G1", NA, x$bee)

x[, treat:= paste(soil, herb, bee, sep = "_")]
x$treat <- ifelse(x$soil == "G1", "G1", x$treat)
setkey(x, treat)

# sep into G1 and treatment gen
g_1 <- x[soil=="G1"]
g_10 <- x[soil!="G1"]

# split by soils
soils <- unique(g_1$soil_plastic)

g_1_l <- g_1[soil_plastic=="L"]
g_1_t <- g_1[soil_plastic=="T"]

# take only g 10 where grown in home soil, remove effect of plasticity
g_10 <- g_10[soil==soil_plastic]

# loop over rep, set g to 10 for ten generations
rep <- c("RA", "RB")
g = 10
rep_list <- vector("list", length = length(rep))

# within in rep, take the trait mean, get the sd and divide by the # of generations
for(l in 1:length(rep)) {
	r_sub <- rep[l]
	print(paste0("loading files for rep: ", r_sub))

	# take only g_1 and g_10 with correct rep
	g_1_rep <- g_1[rep == r_sub]
	
	# get g_1 value depending on soil, remove columns not needed and get trait mean per pop 
	g_1_rep[, c("ind_id","rep", "soil", "herb", "bee", "treat"):=NULL]
	
	# get trait means for g1, depeding on which soil grown in 
	g_1_means <- g_1_rep[, lapply(.SD, mean, na.rm=TRUE), by= c("soil_plastic")]
	traits <- names(g_1_means[,-1])
	g_1_sd <- g_1_rep[, lapply(.SD, sd, na.rm=TRUE), by="soil_plastic"]
	
	# do for g10 
	g_10_rep <- g_10[rep == r_sub]	
	g_10_rep[, c("ind_id","rep", "soil", "soil_plastic",  "herb", "bee"):=NULL]

	# get trait means for g10 split into list and do there
	g_10_means <- g_10_rep[, lapply(.SD, mean, na.rm=TRUE), by=treat]
	g_10_sd <- g_10_rep[, lapply(.SD, sd, na.rm=TRUE), by=treat]
	g10_traits <- colnames(g_10_means)[-1] # get trait names for g10 and remove the treatment col

	mean_list <- split(g_10_means, by=c("treat"))
	sd_list <- split(g_10_sd, by=c("treat"))
	treats <-names(mean_list)
	soils <- str_split_i(treats, "_", 1)
	hald_list <- vector("list", length = length(treats))

	# loop over treatments, easier to read than lapply
	for(i in 1:length(treats)) {
		# get mean per treatment, 
		m_dt <- mean_list[[i]]
		m_dt <- data.table(t(m_dt))
		m_dt <- m_dt[-1,] 
		m_dt <- data.table(cbind(avg = m_dt$V1, trait = g10_traits))
		setkey(m_dt, trait)
		
		# get sd per treatment, create pooled sd
		s_dt <- sd_list[[i]]
		s_dt <- data.table(t(s_dt))
		s_dt <- s_dt[-1,] 
		s_dt <- data.table(cbind(sd = s_dt$V1, trait = g10_traits))
		setkey(s_dt, trait)

		g10 <- m_dt[s_dt]
		setkey(g10, trait)

		# get the g1 values for soil
		s <- soils[i]

		g1_a <- g_1_means[soil_plastic == s]
		g1_s <- g_1_sd[soil_plastic == s]
		g1 <- data.table(cbind(trait = traits, avg = t(g1_a[,-1]), sd = t(g1_s[,-1])))
		setnames(g1, c("V2", "V3"), c("avg", "sd"), skip_absent=TRUE)	
		setkey(g1, trait)

		# subtract from G1, create pooled sd
		g10$dif <- (as.numeric(g10$avg) - as.numeric(g1$avg))
		g10$pooled_sd <- sqrt((as.numeric(g1$sd)**2 + as.numeric(g10$sd)**2)/2)
		
		#g10$pooled_alt <- sqrt(((n1$n-1)*g_1_sd$sd^2 + (n2$V1-1)*as.numeric(s_dt$V1)^2) / (n1$n+n1$n-2)) # cohen alternative, relatively similar tho
		
		# get haldane
		g10$h <- g10[, (dif/pooled_sd)/g]

		# create data table and add to list
		dt <- data.table(cbind(trait = g10$trait, haldane = g10$h, treatment = treats[i]))
		hald_list[[i]] <- dt
	}

	hald <- rbindlist(hald_list)
	hald$rep <- rep[l]
	rep_list[[l]] <- hald
}

full_hald <- rbindlist(rep_list)
full_hald$haldane <- as.numeric(full_hald$haldane) 

fwrite(full_hald, sep = "\t", paste0("haldane_full_rep_soil_sep.tsv"))
