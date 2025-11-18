shared_mods <- function(i) {

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
trait_net <- graph.adjacency(trait_matrix, mode = "undirected")

# add weight, number of shared snps and add layout, not sure if needed
E(trait_net)$weight <- count.multiple(trait_net)
set.seed(123)
l <- layout_with_fr(trait_net) # this layout clusters based on weights, in our case, the number of snps in common

# cluster with three types
print(paste0("time consuming step ~ 10 min: performing clustering for subset ", subset, "- - - - - - - - - - - - - - - - - - - - - - - - -"))
ceb <- cluster_optimal(trait_net)
bet <- cluster_edge_betweenness(trait_net)
eig <- cluster_leading_eigen(trait_net)
print(paste0("finished clustering for subset ", subset, "- - - - - - - - - - - - - - - - - - - - - - - - -"))

# extract and form data frame, get the memebership and see which are shared between the clustering algo
dt <- data.table(cbind(trait = ceb$names, opt_mem = ceb$membership, bet_mem = bet$membership, eig_mem = eig$membership))
n <- choose(nrow(dt),2)

fwrite(dt, sep = "\t", paste0(subset, "_three_modularity_data_frame_sigzone_only.tsv"))

mod_trait <- 1:length(n)
bi_trait <- 1:length(n)
counter = 1
for(k in 1:(nrow(dt)-1)) {
	t1 <- as.character(dt[k,1])
	x1 <- dt[k,2]
	x2 <- dt[k,3]
	x3 <- dt[k,4]
		for(l in (k+1):nrow(dt)) {
			t2 <- as.character(dt[l,1])
			y1 <- dt[l,2]
			y2 <- dt[l,3]
			y3 <- dt[l,4]

		mod_trait[counter] <- paste(t1,t2,sep= "_")
		if (x1==y1 & x2==y2 & x3==y3){bi_trait[counter] <- 1}
		else{bi_trait[counter] <- 0}
		counter = counter + 1
	}
}

# create a data table, 
sub_mod <- data.table(cbind(V1 = mod_trait, bi = bi_trait))
sub_mod$bi <- as.numeric(sub_mod$bi)
sub_mod$trait<- str_split_i(sub_mod$V1, "_", 1)
sub_mod$trait_2<- str_split_i(sub_mod$V1, "_", 2)


# create the matrix
data2 <- data.frame(bi = sub_mod$bi, trait = sub_mod$trait_2, trait_2 = sub_mod$trait) # flip the pop 1 and two 
data3 <- as.data.frame(sub_mod) 
data3<- data3[,-c(1)]
df2 <- rbind(data3, data2)
df3 <- as.data.frame.matrix(xtabs(bi ~ ., df2))
rownames(df3) <- colnames(df3)
df3 <- as.matrix(df3)

write.table(df3)
# get the shared clusters
g1 <- graph.adjacency(df3)
c1 <- clusters(g1)
clust_size <- c1$csize
m <- which(clust_size>1)

# loop over the number of clusters and create the final data set to write out 
print(paste0("final loop for subset ", subset, "- - - - - - - - - - - - - - - - - - - - - - - - -"))
shared_mod <- 1:length(m)
length_mod <- as.character(1:length(m))
final_res <- data.table(cbind(len = length_mod, trait_vector = shared_mod))

for (ik in 1:length(m)) {
	m_check = m[ik]
	mems <- clusters(g1)$mem==m_check
	final_res[ik, 1] <- length(names(c1$membership)[mems])
	final_res[ik, 2] <- paste(names(c1$membership)[mems], collapse = " , ")
}

final_res$len <- as.numeric(final_res$len)
final_res$treat = subset
final_res <- final_res[order(-len)]


return(final_res)

}