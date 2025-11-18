# run simulations of distance matrix to see if we can say two matrixses are dissimilar
if (!require(tidyverse)) install.packages('tidyverse')
if (!require(vroom)) install.packages('vroom')
if (!require(data.table)) install.packages('data.table')
if (!require(igraph)) install.packages('igraph')
if (!require(bipartite)) install.packages('bipartite')
if (!require(scales)) install.packages('scales')
if (!require(svglite)) install.packages('svglite')
#if (!require(bootcluster)) install.packages('bootcluster')
if (!require(geomorph)) install.packages('geomorph')
if (!require(emmeans)) install.packages('emmeans')
if (!require(plot3D)) install.packages('plot3D')
if (!require(fitdistrplus)) install.packages('fitdistrplus')
if (!require(glmmTMB)) install.packages('glmmTMB')
if (!require(DHARMa)) install.packages('DHARMa')
if (!require(car)) install.packages('car')


# to install gg3d
#devtools::install_github("AckerDWM/gg3D")

# loop through treatments
treat <- c("LHB", "LHH", "LNHB", "LNHH", "THB", "THH", "TNHB", "TNHH")

# set wd to place outfile there
file_wd <- "~/mol_eco_2024/2.5.5_gwas/"
out_wd <- "~/mol_eco_2024/2.5.5_gwas/"

dist_list <- vector("list", length = length(treat))
mds_list <- vector("list", length = length(treat))
# pairwise loop to assess network
for(i in 1:(length(treat))) {
	subset = treat[i]
	out_dir <- paste(out_wd, subset, "/", sep = "")
	setwd(out_dir)
	print(paste0("loading files for subset ", subset, "- - - - - - - - - - - - - - - - - - - - - - - - -"))

	# read in treatment file
	pleio_df <- fread(sep = "\t", paste0(subset, "id_trait_pleio_table_sigzone_only.tsv"))
	pleio_df <- unique(pleio_df, by = c('id', "trait"))

	### plot only where pleio tropic, ie 2 or more connections of snp to traits - - - - - - - - - - - - - - - - - - - - - - - - - 
	d_pleio <- pleio_df[pleio >= 2]
	pleio_weights <- unlist(as.vector(abs(d_pleio[, 3])))
	d_pleio <- d_pleio[,1:2]
	d_pleio_mat = table(d_pleio)
	class(d_pleio_mat) <- "matrix" # And we convert it from a table to a matrix

	# transpose to create the trait matrix
	trait_matrix = t(d_pleio_mat) %*% d_pleio_mat
	snp_numb <- diag(trait_matrix)
	diag(trait_matrix) <- 0 # we again set it to 0
	
	trait_matrix <- decostand(trait_matrix, "total")
	dis <-as.matrix(vegdist(trait_matrix,method="bray"))
	dist_list[[i]] <- dis

	mds <- cmdscale(dis, k = 3)
	mat <- as.matrix(cbind(x = mds[,1], y = mds[,2], z = mds[,3]))
	mds_list[[i]] <- mat
}

# plot mds
#plot(x,y)
#text(x,y,cex= 0.7, labels=row.names)


# move distances into array
column.names <- colnames(dis)
row.names <- rownames(dis)
matrix.names <- paste("dis_", treat, sep = "")

# Take these vectors as input to the array.
res_arr <- array(unlist(dist_list), dim = c(nrow(dis),nrow(dis), length(dist_list)),
		dimnames = list(row.names,column.names,matrix.names))


# convert to geomorph data frame to add the factor levels
data <- geomorph.data.frame(mat = res_arr, soil = str_sub(treat,1,1), herb = str_sub(treat,2,-2), bee = str_sub(treat,-1,-1)) 


Y.gpa <- gpagen(res_arr)


# test for differences with mds instead - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# move them into array
column_names_mat <- colnames(mat)
row_names_mat <- rownames(dis)
matrix_names_mds <- paste("mds_", treat, sep = "")

# Take these vectors as input to the array.
mds_arr <- array(unlist(mds_list), dim = c(nrow(mds),ncol(mds), length(mds_list)),
		dimnames = list(row_names_mat,column_names_mat,matrix_names_mds))


# convert to geomorph data frame to add the factor levels
mds_data <- geomorph.data.frame(mat = mds_arr, soil = str_sub(treat,1,1), herb = str_sub(treat,2,-2), bee = str_sub(treat,-1,-1)) 

Y.gpa <- gpagen(mds_arr)
y<-two.d.array(Y.gpa$coords)


### run anova to see if networks differ by shape, then by size, does the interaction match the genetic or the phenotypic divergence
lm_shape_pheno <- procD.lm(Y.gpa$coords ~ soil*bee, data = mds_data ,iter=500000, int.first = TRUE) 
summary(lm_shape_pheno)

#          Df      SS      MS     Rsq      F        Z Pr(>F)
#soil       1 0.45145 0.45145 0.17327 1.2458  0.94133 0.1620
#bee        1 0.27768 0.27768 0.10658 0.7663 -0.56995 0.7546
#soil:bee   1 0.42677 0.42677 0.16380 1.1777  0.30231 0.3434
#Residuals  4 1.44957 0.36239 0.55635                       
#Total      7 2.60548                                                                           


lm_shape_gen <- procD.lm(Y.gpa$coords ~ herb*bee, data = mds_data ,iter=500000, int.first = TRUE) 
summary(lm_shape_gen)
                             
#          Df      SS      MS     Rsq      F        Z Pr(>F)
#herb       1 0.40769 0.40769 0.15647 1.0007  0.14036 0.4637
#bee        1 0.27768 0.27768 0.10658 0.6816 -0.84621 0.8618
#herb:bee   1 0.29051 0.29051 0.11150 0.7131 -0.50587 0.7248
#Residuals  4 1.62959 0.40740 0.62545                       
#Total      7 2.60548      


### overall no differences between shape, but what about size
# makes sense to convert size to data frame to just run linear model
size_dt <- data.table(cbind(size = Y.gpa$Csize, soil = mds_data$soil, herb = mds_data$herb, bee = mds_data$bee))
setwd(out_wd)
fwrite(size_dt, sep = "\t", paste0("network_size_treatments_sigzone_only.tsv")) 

descdist(as.numeric(size_dt$size), discrete = FALSE) 
# model
lm_size_pheno <- glmmTMB(1/(as.numeric(size_dt$size)) ~ soil*bee, data = size_dt, family = beta_family())
Anova(lm_size_pheno)
#Response: 1/size_dt$size
#           Chisq Df Pr(>Chisq)    
#soil      4.4024  1  0.0358890 *  
#bee       5.6106  1  0.0178517 *  
#soil:bee 11.9393  1  0.0005496 ***


summary(lm_size_pheno)
#            Estimate Std. Error z value Pr(>|z|)    
#(Intercept)  0.006882   0.027687   0.249    0.804    
#soilT        0.056342   0.039165   1.439    0.150    
#beeH         0.043828   0.039161   1.119    0.263    
#soilT:beeH  -0.230590   0.055411  -4.161 3.16e-05 ***
      

lm_size_gen <- glmmTMB(1/(as.numeric(size_dt$size)) ~ herb*bee, data = size_dt, family = beta_family()) 
summary(lm_size_gen)

#             Estimate Std. Error z value Pr(>|z|)
#(Intercept) -0.035428   0.045238  -0.783    0.434
#erbNH       0.064924   0.063974   1.015    0.310
#eeH         0.001365   0.063975   0.021    0.983
#herbNH:beeH -0.144234   0.090507  -1.594    0.111

Anova(lm_size_gen)
#Response: 1/(as.numeric(size_dt$size))
#          Chisq Df Pr(>Chisq)
#herb     0.0249  1     0.8746
#bee      2.4409  1     0.1182
#herb:bee 2.5396  1     0.1110  

## combine to big table
a <- summary(lm_shape_pheno)
a_1 <- data.table(a$table)

a_g <- summary(lm_shape_gen)
a_g_1 <- data.table(a_g$table)

s_1 <- data.table(anova(lm_size_pheno))

s_g_1 <- data.table(anova(lm_size_gen))


table_list <- list(a_1, a_g_1)
table_list_s <- list(s_1, s_g_1)

t_1 <- rbindlist(table_list) 
t_2 <- rbindlist(table_list_s) 

#setwd(out_wd)
#fwrite(t_1, sep = ",", paste0("anova_shape.csv")) 
fwrite(t_2, sep = ",", paste0("anova_size.csv")) 

size_dt[, round(range(as.numeric(size)), 2), by = c("soil", "bee")]
size_dt[, round(mean(as.numeric(size)), 2), by = c("soil", "bee")]
size_dt[, round(sd(as.numeric(size)), 2), by = c("soil", "bee")]

# so there seems to be an interaction for size, tuff soil with hand tends to be bigger, thus wider and more spread apart, less connectivity
size_em <- emmeans(lm_size_pheno, ~  soil*bee, type = "response"); size_em
size_em <- data.table(as.data.frame(size_em))
fwrite(size_em, sep = "\t", paste0("network_size_soil_bee.tsv")) 
