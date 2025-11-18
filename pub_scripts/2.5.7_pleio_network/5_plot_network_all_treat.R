# create 3d plot, sigzone only
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
if (!require(viridis)) install.packages('viridis')

# to install gg3d
#devtools::install_github("AckerDWM/gg3D")

# loop through treatments
treat <- c("LHB", "LHH", "LNHB", "LNHH", "THB", "THH", "TNHB", "TNHH")

# set wd to place outfile there
out_wd <- "~/mol_eco_2024/2.5.5_gwas/"

dist_list <- vector("list", length = length(treat))
mds_list <- vector("list", length = length(treat))

# read in size file, may be best to structure as largest to smallest networks to see distinction
setwd(out_wd)
size_dt <- fread(sep = "\t", paste0("network_size_treatments_sigzone_only.tsv"))

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
	
	dis <-as.matrix(vegdist(trait_matrix,method="bray"))
	dist_list[[i]] <- dis

	mds <- cmdscale(dis, k = 3)
	mat <- as.matrix(cbind(x = mds[,1], y = mds[,2], z = mds[,3]))
	mds_list[[i]] <- mat
}



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


# create color scale for trait class, and for easier membership loading - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
setwd(out_wd)
trait_df <- fread(sep = "\t", paste0("~/mol_eco_2024/2.5.5_gwas/trait_class.tsv")) 
trait_df$trait_class <- as.factor(trait_df$trait_class)

col_max <- length(unique(trait_df$trait_class))
col <- data.table(color = viridis_pal(alpha = 1, option = "turbo")(col_max))
col = cbind(col, trait_class = unique(trait_df$trait_class))
setkey(col, "trait_class")
setkey(trait_df, "trait_class")

trait_df <- col[trait_df] # merge color then set up to merge with mds
setkey(trait_df, "trait")

# load in treatment modules membership and plot based on that - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
mem_total <- fread(sep = "\t", paste0("treat_mod_mem_sigzone_only.tsv"))
mem_total <- data.table(t(mem_total))
colnames(mem_total) <- treat

# get max number of modules and create color based on that
col_max_mod <- max(mem_total, na.rm = T)
col_mod <- data.table(color = viridis_pal(alpha = 1, option = "turbo")(col_max_mod))

col_mod = cbind(color_modules = col_mod, mem = 1:col_max_mod)
setkey(col_mod, "mem")

# get plotting data frame for all treatments and add to list for easier opperations
treat_list <- vector("list", length = length(treat))

for (i in (1:length(mds_list))) {
	dt <- data.table(mds_list[[i]])
	
	# merge with trait name and for color
	dt$trait <- row_names_mat
	setkey(dt, "trait")
	dt <- trait_df[dt]

	# merge with mod membership
	dt$mem <- mem_total[,..i]
	setkey(dt, "mem")
	dt <- col_mod[dt]

	# add to list and name based on treat
	treat_list[[i]] <- assign(treat[i], dt)
	
}


# get min and max so on same scale
mx <-min(unlist(lapply(treat_list, function(dt) {min(dt$x)})))
my <- min(unlist(lapply(treat_list, function(dt) {min(dt$y)})))
mz <-min(unlist(lapply(treat_list, function(dt) {min(dt$z)})))

max_x <- max(unlist(lapply(treat_list, function(dt) {max(dt$x)})))
max_y <-max(unlist(lapply(treat_list, function(dt) {max(dt$y)})))
max_z <-max(unlist(lapply(treat_list, function(dt) {max(dt$z)})))

# create rank to plot, smallest to biggest
p_r <- rank(size_dt$size)

# plot big boy
resolution <- 144
svglite(paste0("all_treat_three_d_mod_traits_sigzone_only.svg", sep=""), width = 1080/resolution, height = 648/resolution)
#tiff(paste0("three_d_mod_traits.tiff", sep=""), res=1200, width=25, height=15, unit="cm",compression="lzw") 

# set plotting window
par(mfrow = c(2, length(treat)/2),
	mai = c(0.0, 0.3, 0.3, 0.2)) # to reduce space between plts

for (i in (1:length(treat_list))) {
	# plot based on network size, smallest to largest
	j <- which(p_r==i)
	dt <- treat_list[[j]] 

	scatter3D(dt$x, dt$y, dt$z, pch = 19, 
		xlim=c(mx,max_x), ylim=c(my,max_y), zlim=c(mz,max_z), type = "h", phi = 10,
		colkey = FALSE, colvar = as.integer(dt$mem), col = col_mod$color_modules, NAcol = "#9E9E9E",
		main = paste0(treat[j], " - Density:" , round(1/size_dt$size[j], 2)))
	#text3D(dt$x, dt$y, dt$z,  labels = dt$trait, adj = -0.125,
	        #add = TRUE, colkey = FALSE, cex = 0.7)

}

dev.off()