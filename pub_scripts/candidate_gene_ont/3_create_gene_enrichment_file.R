# perform gene enrichment
if (!require(data.table)) install.packages('data.table')
if (!require(GenomicRanges)) install.packages('GenomicRanges')
if (!require(tidyverse)) install.packages('tidyverse')
if (!require(BiocManager)) install.packages('BiocManager')

# install cluster
BiocManager::install("clusterProfiler")

# require cluster
if (!require(clusterProfiler)) install.packages('clusterProfiler')

# read in gene file
setwd("~/mol_eco_2024/2.5.6_gene_ont")
full_gene <- fread(sep = "\t", paste0("gene_annotation_treat_loop_shared_clear_cmh.tsv"))

full_gene <- fread(sep = ",", paste0("gene_annotation_treat_clear_cmh_gwas_supp_table.csv"))

# read in gene file 
bra_go <- read.delim("~/mol_eco_2024/2.5.6_gene_ont/transcripts_go_exp6879.txt", header = TRUE)

# make name and gene links
term2gene <- bra_go[,c(3,2)]
all_genes <- unique(term2gene$transcript_id)
term2name <- unique(bra_go[,c(3,6)])

# define treatments 
treat <- c("LHB", "LHH", "LNHB", "LNHH", "THB", "THH", "TNHB", "TNHH")
treat_list <- vector("list", length = length(treat))

# loop over treatment
for(i in 1:length(treat)) {
  dt <- as.data.table(full_gene[full_gene$treatment == treat[i],])
  setkey(dt, ID)

  # read in process name file
  bra_sub_list <-dt$gene

  # test for enrichment 
  go_out <- enricher(
    bra_sub_list,
    pvalueCutoff = 0.99,
    pAdjustMethod = "BH",
    all_genes,
    TERM2GENE = term2gene,
    TERM2NAME = term2name,
    qvalueCutoff = 0.9
  )

  # save as data frame, then data table and then save to list
  go_out_df <- as.data.frame(go_out)
  go_out_df$treat <- treat[i]
  treat_list[[i]] <- as.data.table(go_out_df)
}

# bind list to data table 
total_gene_ont_dt <- rbindlist(treat_list)

# write out file
fwrite(total_gene_ont_dt, sep = "\t", paste0("~/20241015_total_gene_ontology.tsv"))