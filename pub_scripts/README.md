# Scripts for: Soil and Herbivory Alters Genomic Adaptation to Bumblebee Pollination 
We used genomic and phenotypic data to better understand how ecological factors affect genomic divergence between three factors: soil, herbivory, and pollination. 

The genomic sequences analyzed in this article are stored and accessible through the National Center for Biotechnology Information (NCBI) under https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1105729. 

R scripts and metadata are stored on dryad: https://doi.org/10.5061/dryad.280gb5n0x

# Description of the data and file structure
These scripts will recreate the analyses shown in the publication. The heading numbers are based on the material methods section. Within in the folders, scripts are structured numerically. Thus, 1_ is first, 2_ is second and so on. You may need to change the working directory to rerun the scripts. 

The archive.tar.gz contains a combination of raw and metadata. Please see PLINK, EMMAX and the manuscript or supplement for further documentation.

## Treatment info 
Soil: L = limestone soil
T = tuff soil
Herbivory: H = herbivory
NH = no-herbivory
Bee / Pollination: B = bee pollination
H = hand pollination

# Headers:
# 2.4_pheno_evo - Creates phenotypic trait evolution analyses
## file: emmax_pheno_data.tsv
purpose: to run phenotypic genomic associations, 
comments on missing data: volatile data that is NA is because there was not enough volatile to register after the threshold or 
because there were not enough flowers to measure volatiles, 
other data that is missing is because plants did not produce traits on the day of measurement
id column description: ind_id = plant id, 
morphological column descriptions: heightd22 and heightd34 = plant height at day 22 and 34 after sewing in cm, 
leaf-size = (leaf widht x leaf length)/2 in cm,
flowering time = days after sewing for first flower, branches lenght = average length of branches in cm, nbbranches = number of branches,
flowersperdaysincefloweringtime = number of flowers per day after starting flowering, flowerdiamter = flower diameter in cm,
petal-length and petal-width = average measured in cm, sepale-length and style and stamen = measured in cm
floral volatiles column descriptions: all volatiles are measured in pg/flower/hour and natural log transformed

## file: haldane_full_rep_soil_sep.tsv
purpose: to run model assessing haldane
column descriptions: trait = the trait for the haldane being calculated, haldane = the haldane calculated generation 10 - generation 1,
treatment = which treatment, rep = for which replicate

## file: n_snps_clear_cmh_gwas_haldane_full_rep_soil_sep.ts
purpose: to run model assessing assocation between haldane and number of selected snps
column descriptions: trait = the trait for the haldane being calculated, avg_haldane = average haldane per treatment,
treatment = which treatment, n_snps = the number of selected snps indentified in selection scans

## file: shared_cmh_af_all_treatments_traits.tsv
purpose: shows the allele frequency, pvalue for snps, and associations with traits 
column descriptions: chr = chromosome, pos = position on chromosome, cmh_fdr = fdr corrected pvalue of cmh tests,
G10a_af = the allele frequency of the last generation for the allele 1, G1a_af = the allele frequency of the first generation for the allele 1, 
id = chromosome and pos combined or snp id, trait = phenotypic trait, 
treatment = which treatment, af_diff = diffence in allele frequency between generation 10 and 1,
trans_af_diff = transformed af_diff with arcsine square root

# 2.5.1_fst_pair - Fst files to run pairwise fst analyses
## file: example : RA_G1_RB_G1_pop.windowed.weir.fst
File is pairwise fst between two populations. RA and RB are replicates (A or B). After the underscore is the population. 
In this example, it is between both generation 1 populations. 
column descriptions: CHROM = chromosome, BIN_START = start of 5000 bp windows, BIN_END = end of 5000 bp windows,
N_VARIANTS = number of snps within window, WEIGHTED_FST = mean fst weighted by number of variants,
MEAN_FST = mean fst calculated for window

## file: example : RA_LLNHB_rep_pop
rep_pop files contain which plants are in the population and used for fst analyses

## file: example : RA_fst_full_data.tsv
It is the aggregation of fst results for that replicate (A or B)

## file: total_rep_sep_fst_full_data.tsv
It is the aggregation of the two full_data files above 

## file: fst_between_g1_and_all_treat.tsv
It has fst between g1 and that treatment for each replicate, and averaged over replicate (m column)
It is the aggregation of fst results for that replicate (A or B)

# 2.5.2 Create Fst model, and 1000 permutations
## file: fst_rrpp_perm_boxcox_g10_only_trans.tsv
contains model results (F values) for each fixed effect and random effect

# 2.5.3 Runs selection scans for CMH tests and CLEAR
## file: example : dr_RA_G1.frq
contains allele frequency for the population. Here is it is replicate RA (A) and population generation 1.
column descriptions: CHROM = chromosome, POS = bp position on chromosome, 
N_ALLELES = number of alleles, N_CHR = number of chromosomes (copies essentially), 
{ALLELE:FREQ} = The base and allele frequency of both alleles 1 and 2

## file: clear_and_cmh_all_treat_pruned.tsv
contains cmh_fdr and CLEAR s to run model and assess signifigance for each SNP for each treatment
column descriptions: fdr = fdr adjusted pvalue, s = selection coefficient from CLEAR, 
sig = if both fdr and s are above threshold used to assess significance

## file: model_both_clear_cmh_num_snps_soil_sep.tsv
model results from above file with mean model estimate, lower and upper confidence intervals and actual number of snps. Treatment factors are listed
column descriptions: mean = mean model estimate, se = standard error, 
full_group = numbers indicating pairwise signifigance between all treatments,
soil_group = numbers indicating pairwise signifigance between all treatments in same soil,
lcl and ucl = 95% confidence intervals

## file: shared_snp_exact_pos_mat.csv
matrix indicating the number of shared selected snps between treatments

## file: example : subsetG1_A.prune.in
contains the snps which after prunning are used for testing, done to account for linkage disequilibrium

## files: example: subset_LHBcmh_fdr_sig_updated_ne_pruned_prior.tsv
Contains the allele frequencies and pvalues for that specific treatment for the cmh test
G10 is generation 10
G10a is generation 10 replicate A
af = allele frequency, base = nucleic base

## files: avg_ne.csv
contains effective population size to run clear and cmh tests
column descriptions: V1 = effective population size, rounded_ne = ne rounded up

# 2.5.4 Creates genomic PCAs
We ran two PCAs 1) full data set and 2) on selected snps only
.eigenval suffix contains eigenvals for the PCA.  
.eigenvec contains the PCA values for each individual plant all the way to PC10

# 2.5.5 Runs Genome wide assocation studies (GWAS)
.bed, .tfam, .tped, are all plink generated files to get genotype data, please see plink for more information
suffix .aIBS.kinf files are kinship matrixes generated in emmax. Please see emmax for better documentation.

The data subdirectories (LHB, LHH) contain the gwas files for each treatment and for each phenotypic trait. 
the files with the suffix: emmax.reml contain the emmax maximum likelihood estimation, please see emmax for better documentation.

## files: example: LHB_Sepale-length_emmax_localscore_full.tsv
contains the phenotypic association for that chromosome and position and if it is a significant zone
column descriptions: beta = the effect size of the estimated assocation,
lindley and L and score, and cor and th = used in local score to assess signficant zones,
mahattan = is the bp when plotting, col = color when plotting, sz_filter = is the snp in a significant zone 

## files: example: LHB_Petal-width_emmax_fdr_sig.tsv
contains the phenotypic association for that chromosome and position and only if it is a significant zone

## files: example: LHB_Stamen_emmax_only_sigzones_1kb_range.tsv
contains the significant zone for the phenotypic association plus/minus 1kb 
column descriptions: beg and end = start and end of significant zone, peak = how tall is the peak in the zone,
QTL_length = length of zone, minus_kb and plus_kb = start and end of zone plus minus 1000,
zone_id = chr and end pasted to get a zone name

## files: example: LHB_Style_emmax_sig_snps_sigzones.tsv
contains all the signicaint snps within a significant zone

## files: example: LHB_pleio_table_sigzone_only.tsv
contains the number of snps exhibiting pleiotropy for that treatment and level of pleiotropy

## files: example: LHB_shared_modules_sigzone_only.tsv
contains the modules for that treatment
column descriptions: len = number of traits in module, trait_vector = which traits are in the module, treat = treatment

## files: sup_mod_list_sigzone_only.csv
aggregation of file directly 1 above

## files: example: LHB_three_modularity_data_frame_sigzone_only.tsv
contains the modules for that treatment assessed with three methods, please see manuscript for further details
Numbers tell which module number the trait is for that membership algorithim

## files: example: LHBid_trait_pleio_table_sigzone_only.tsv
contains the zones which exhbit pleiotropy and for which traits

## files: example: _phen_emmax.tsv
contains the phenotypes for that specific treatment, see above,

## files: emmax_phenotypic_data_thomas.csv
contains the phenotypes for that specific treatment, see above,

## files: trait_class.tsv
trait_class = whether trait is morphological or floral volatile,

## files: emmax_estimates_all_treat_trait.tsv
has the variance calculated by emmax for each treatment and trait, please see emmax for better documentation
column descriptions: trait_class = whether trait is morphological or floral volatile, 
va = additive variance, ve = environmental variance, va_ve = ratio of va/ve, pseudo_h = psuedo heritability

## files: network_size_treatments_sigzone_only.tsv
contains the network size generated from pleiotropy

## files: pleio_index_sigzone_only.csv
contains the pleiotropy index for each treatment, please see manuscript for calculations

## files: shared_snps_between_traits_sigzone_only.tsv
for each treatment and trait pair, it contains the number of shared snps associated with both traits

## files: treat_mod_mem_sigzone_only.tsv
contains matrix of which traits are in which modules

# 2.5.6 Creates pleiotropy analyses based on GWAS. 
# 2.5.7 Creates pleiotropic trait networks for each treatment
# 2.5.8 Creates modules within pleiotropic trait networks

# 2.5.6_gene_ont
## files: GCF_000309985.2_CAAS_Brap_v3.01_genomic_converted.txt
contains which candidate genes were annotated to the region of the genome
column descriptions: type = whether it is a gene, exon or mRNA among others, gene = gene name,
product and biotype = function of the gene

## files: transcripts_go_exp6879.txt
contains go terms of genes
column descriptions: go = goterm, 

## files: shared_clear_cmh_gwas_all_treatments_traits.tsv
contains whether the genomic region is significant and associates with phenotypic traits

# Sharing/Access information
The raw phenotypic and genomic reads data used in this publication is published in a previous study. 
Genomic metadata from this publication is generated for this publication. 