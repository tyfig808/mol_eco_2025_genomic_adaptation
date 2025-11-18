# gene annotation with gene retriever
if (!require(data.table)) install.packages('data.table')
if (!require(GenomicRanges)) install.packages('GenomicRanges')

setwd("~/mol_eco_2024/2.5.5_gwas/")

# loop over treatments 
treat <- c("LHB", "LHH", "LNHB", "LNHH", "THB", "THH", "TNHB", "TNHH")

# load gene file
genes <- fread(sep = "\t", paste0("~/mol_eco_2024/2.5.6_gene_ont/GCF_000309985.2_CAAS_Brap_v3.01_genomic_converted.txt"))
df_genes <- genes[type == "exon"]

# convert to numerical chrome
x1 <- c(1,2,3,4,5,6,7,8,9,10)
x2 <- c("A01","A02","A03","A04","A05","A06","A07","A08","A09","A10")
look <- data.table(x1,x2)
setkey(look, x2)

# merge to get numerical chrome
setkey(df_genes, chromosome)
df_genes <- df_genes[look]
setnames(df_genes, c("x1"), c("chr"), skip_absent=TRUE)	

setkey(df_genes, gene)

# read in file you want to know SNPs for, this is where snps are in cmh and gwas
freq <- fread(sep = "\t", paste0("shared_clear_cmh_gwas_all_treatments_traits.tsv")) 
freq <- na.omit(freq, "chr") 


treat_list <- vector("list", length = length(treat))

for(i in 1:length(treat)) {
	dt <- freq[treatment == treat[i]]
	setkey(dt, id)
	# use genomic ranges to get subset of hits
	win = 1000 # Windows 1kb up and downstream

	gr1 <- with(dt, GRanges(chr, IRanges(start=pos, end=pos, names=id), SNP=id))
	gr2 <- with(df_genes, GRanges(chr, IRanges(start=start-win, end=end+win, names=gene), ATG=gene))
	ranges <- subsetByOverlaps(gr2,gr1)
	hits <- findOverlaps(gr2, gr1)

	# take id from our subset, called rsid, and add close gene in ATG col
	rsid <- CharacterList(split(names(gr1)[subjectHits(hits)], queryHits(hits))) 
	mcols(ranges) <- DataFrame(mcols(ranges), rsid)
	ranges

	df <- data.table(data.frame(mcols(ranges)))
    setnames(df, c("ATG", "rsid"), c("gene", "ID"), skip_absent=TRUE)   

    # merge to get gene product and biotype
    setkey(df, gene)
    df <- df[df_genes, nomatch = NULL]
    df<-df[!duplicated(gene)]

    # list of ID so loop over list and take first element as ID
    df$ID <- sapply(df$ID, "[[", 1)
    setkey(df, ID)

    # merge and drop the dupes of id and trait
    df <- dt[df, nomatch = NULL]
    df<-unique(df, by=c("id", "trait"))

	df$treatment <- treat[i]

	treat_list[[i]] <- df

}

# create full file, remvoe extra col and write out
full_gene <- rbindlist(treat_list)
full_gene[, c("G10a_af", "G1a_af", "type", "i.chr"):=NULL] 



fwrite(full_gene, sep = "\t", paste0("gene_annotation_treat_loop_shared_clear_cmh_gwas.tsv"))

# write out for excel in sup
supp_gene <- full_gene
supp_gene[, c("af_diff", "trans_af_diff", "chromosome", "biotype", "start", "end"):=NULL]

fwrite(supp_gene, sep = ",", paste0("gene_annotation_treat_clear_cmh_gwas_supp_table.csv"))