module load anaconda3
eval "$(conda shell.bash hook)"
source activate vcftools

wd=~/mol_eco_2024

# from full filter data set, create fdr subset
vcftools --vcf final_filter_no_indel_two_allele.recode.vcf --bed ~/mol_eco_2024/
clear_cmh_all_treat_sig_only_pruned.bed --out clear_cmh_sig_no_indel --recode --keep-INFO-all

# create pca for clear and cmh subset
plink2 --vcf clear_cmh_sig_no_indel.recode.vcf --out clear_cmh_sig_no_indel --make-pgen --allow-extra-chr
plink2 --pfile clear_cmh_sig_no_indel --pca --allow-extra-chr --out pca_clear_cmh_sig_no_indel

# create pca for full
plink2 --vcf final_filter_no_indel_two_allele.recode.vcf --out final_total_snp_no_indel --make-pgen --allow-extra-chr
plink2 --pfile final_total_snp_no_indel --pca --allow-extra-chr --out pca_final_total_snp_no_indel