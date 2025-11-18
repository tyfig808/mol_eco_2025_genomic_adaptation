### as calc the sigzones can freeze up if there are peaks in the beginning, this accounts for this and is more robust, ie will not freeze
### gives the same results when the early peak is not there !!!
### only a little slower so can use this from the start or as needed

setwd("~/mol_eco_2024/2.5.5_gwas")

if (!require(tidyverse)) install.packages('tidyverse')
if (!require(vroom)) install.packages('vroom')
if (!require(data.table)) install.packages('data.table')
if (!require(qqman)) install.packages('qqman')

### run pheno file to get the number of cols
phen_files <- list.files(pattern = "_phen_emmax.tsv", recursive = TRUE, full.names = TRUE)

## for tfam files, need to create emmax pheno file, although the emmax wants the tped prefix, same number of characters so it works
t_files <- list.files(pattern = "_G1.numchr.tped", recursive = TRUE, full.names = TRUE)

f_names <- str_split_i(phen_files, "_", 1)
f_names <- str_sub(f_names,3,6)

# for tfam files, emmax just wants the prefix
t_names <- str_split_i(t_files, "/", 2)
t_names<- str_sub(t_names,1,-6)

# use this to run with array in bash at same time
args <- commandArgs(trailingOnly=T)
j <- eval(parse(text=args[1]))
#j = 3

subset <- f_names[j]
print(paste0("loading files for subset ", subset))

phen <-  vroom(phen_files[j]) # think we need the head here

# set wd to place outfile there
file_wd <- "~/mol_eco_2024/2.5.5_gwas"
out_wd <- "~/mol_eco_2024/2.5.5_gwas"
out_dir <- paste(out_wd, subset, "/", sep = "")
setwd(out_dir)

# read in emmax files, loop over
ps_files <- list.files(pattern = "emmax.ps", recursive = TRUE, full.names = TRUE)
ps_names <- str_split_i(ps_files, "[.]", 2)
ps_names <- unlist(regmatches(ps_names,gregexpr("(?<=_).*",ps_names,perl=TRUE))) # extract after the first _

tped <-  fread(paste("/", t_files[j], sep = ""), select = c(1,4))
setnames(tped, c("V1", "V4"), c("CHROM", "POS"), skip_absent=TRUE)	

#Import functions.
source('~/mol_eco_2024/scorelocalfunctions.R')
source('~/mol_eco_2024/peak_beg_chr.R') # load in tyler correction for peak at beginning

for(i in 1:length(ps_files)) {
	#i <- 14	
	trait <- ps_names[i]
	print(paste0("loading files for subset ", subset, " and for trait: ", trait))


	# load phenotypic data and add the chr and pos each time, then write it, emmax then loads it in
	res <- fread(ps_files[i], select = c(2:4))
	setnames(res, c("V2", "V3", "V4"), c("beta", "se", "pvalue"), skip_absent=TRUE)
	res <- cbind(tped, res)
	res <- na.omit(res) 

	# create new cols for fdr and log pvalue 
	res[,log_pval:=-log10(res$pvalue)]
	res[,fdr:=p.adjust(res$pvalue, method = "fdr")]

	# order
	keycol <-c("CHROM","POS")
	res<- setorderv(res, keycol)


	setkey(res, CHROM)
	Nchr=length(res[,unique(CHROM)])

	res[pvalue ==0] # check if there any 0, should not be but if there are then do this code below
	#min(res$pvalue[-which(res$pvalue ==0)])

	chrInfo=res[,.(L=.N,cor=autocor(pvalue)),CHROM]
	setkey(chrInfo,CHROM)
	data.table(chr=res[,unique(CHROM),], S=cumsum(c(0,chrInfo$L[-Nchr])))

	# The score mean must be negative, ksi must be chosen between mean and max of -log10(p-value)
	# xi should be between 2 and 3 for gwas, make sure that mean score is negative
	xi=2
	res[,score:= res$log_pval-xi]
	mean(res$score)
	res[,lindley:=lindley(score),CHROM]

	mean(res$log_pval);max(res$log_pval)

	#hist(mydata$score)
	max(res$lindley)

	# Compute significance threshold for each chromosome
	## Uniform distribution of p-values
	chrInfo[,th:=thresUnif(L, cor, xi),CHROM]
	(res=res[chrInfo])

	# loop over chrome for sigzone
	ch <- unique(res$CHROM)
	sig_list <- vector("list", length = length(ch))

	### still have problems with function looping over more than chrome so we loop over chrome just for sigzone
	for (k in 1:length(ch)){  
		# seperate by chrome, add chrome col, rearrange order and add to list
		ch1 = res[res$CHROM == ch[k]]
		sigZones=sig_sl_chr_8(ch1$lindley, as.numeric(ch1$POS), unique(ch1$th))

		sigZones$CHROM <- ch[k]
		setcolorder(sigZones, c(4,1,2,3)) # move column up
			
		sig_list[[k]] <- sigZones
		print(paste("chromesome: ", k, " completed", sep = ""))
	}

	# bind back and do operatations
	sigZones <- rbindlist(sig_list) 

	sigZones$QTL_length <- sigZones$end - sigZones$beg
	sigZones <- filter(sigZones,peak>0)
	sigZonesSNPDens = data.frame()
	all_SNP_signif_log4 = vector()

	# round the cols
	res$signif = round(res$lindley-res$th,3)
	res$th = round(res$th,3)
	res$cor = round(res$cor,3)
	res$log_pval = round(res$log_pval,6)

	for (zone in c(1:nrow(sigZones))) {
	 CHR = sigZones$chr[zone]
	 BEG = sigZones$beg[zone]
	 END = sigZones$end[zone]
	 SNP_number_log4 = nrow(ch1[ch1$CHROM==CHR & ch1$POS>=BEG & ch1$POS<=END & ch1$log_pval>4])
	 SNP_signif_log4 = ch1$rs[ch1$CHROM==CHR & ch1$POS>=BEG & ch1$POS<=END & ch1$log_pval>4]
	 SNP_number_log5 = nrow(ch1[ch1$CHROM==CHR & ch1$POS>=BEG & ch1$POS<=END & ch1$log_pval>5])
	 SNP_number_log6 = nrow(ch1[ch1$CHROM==CHR & ch1$POS>=BEG & ch1$POS<=END & ch1$log_pval>6])
	 subsetzone = cbind(sigZones[zone,],SNP_number_log4,SNP_number_log5,SNP_number_log6)
	 sigZonesSNPDens = rbind(sigZonesSNPDens,subsetzone)
	 all_SNP_signif_log4 = append(all_SNP_signif_log4,SNP_signif_log4)
	}
	

	# order sigzones
	keycol <-c("CHROM","beg")
	sigZones<- setorderv(sigZones, keycol)	

	
		#write sigzone to table
	fwrite(sigZones, sep = "\t", paste0(subset, "_", trait, "_emmax_sigzones.tsv")) #cmh local score with sig zones here 

	### create manhattan col for plotting and save
	cmh <- as.data.frame(res)
	cmh$manhattan=0
	colnames(cmh)[c(1,2)] =c("chr", "pos")
	for (c in 1:length(unique(cmh$chr))){
	  assign(paste("cmh",c,sep=""),cmh[cmh[,1]==unique(cmh$chr)[c],c(1:ncol(cmh))]) # create data frame
	}

	cmh1$manhattan=cmh1$pos
	cmh2$manhattan=cmh2$pos+max(cmh1$manhattan)
	cmh3$manhattan=cmh3$pos+max(cmh2$manhattan)
	cmh4$manhattan=cmh4$pos+max(cmh3$manhattan)
	cmh5$manhattan=cmh5$pos+max(cmh4$manhattan)
	cmh6$manhattan=cmh6$pos+max(cmh5$manhattan)
	cmh7$manhattan=cmh7$pos+max(cmh6$manhattan)
	cmh8$manhattan=cmh8$pos+max(cmh7$manhattan)
	cmh9$manhattan=cmh9$pos+max(cmh8$manhattan)
	cmh10$manhattan=cmh10$pos+max(cmh9$manhattan)

	c_list=list(cmh1 = cmh1,cmh2 = cmh2,cmh3 = cmh3,cmh4 = cmh4,cmh5 = cmh5,
	  cmh6 = cmh6,cmh7= cmh7,cmh8 = cmh8,cmh9 = cmh9,cmh10 = cmh10)

	res <- as.data.table(rbindlist(c_list))


	### create colors and min max 
	res$col <- ""

	odd <- c(1,3,5,7,9)  
	even <- c(2,4,6,8,10)  

	res$col <- ifelse(res$chr %in% odd, 
	                       "grey20", 
	                       ifelse(res$chr %in% even, "grey40", NA))

	deb_emmax=min(res$manhattan)
	fin_emmax=max(res$manhattan)

	### extract snps from sigzone and fdr
	keycol <-c("CHROM","beg")
	sigZones<- setorderv(sigZones, keycol)
	#setDT(res)[, sz_filter := FALSE]

	res <- res[sigZones, sz_filter := chr == CHROM, # first match same chrome
	       on = .(pos >= beg, pos <= end)] # then check if bigger than begining and smaller then end

	# interestingly, filtering with the sz filter and snp4 give different snps, should go with the snp4 because it is the way they designed it
	#snps_sigzones <- res[res$sz_filter == TRUE] 
	#SNP_signif_log4<- res[signif>0]
	SNP_signif_log4<- res[signif>0 & log_pval>4]

	fdr_sub <- res[fdr<0.05] # subset with only sig fdr 
	fwrite(fdr_sub, sep = "\t", paste0(subset, "_", trait, "_emmax_fdr_sig.tsv")) # fdr correction and snps sig after

	# write the full data set 
	fwrite(res, sep = "\t", paste0(subset, "_", trait, "_emmax_localscore_full.tsv")) # full data with distances

	# write the sigzones
	fwrite(SNP_signif_log4, sep = "\t", paste0(subset, "_", trait, "_emmax_sig_snps_sigzones.tsv")) #cmh local score with sig zones here 

	### plot now
	tiff(paste0(subset, "_", trait, "_emmax_snps_sigZones.tiff", sep=""), res=600, width=25, height=15, unit="cm",compression="lzw") 

	m=matrix(1:3,3,1)
	layout(m,c(5),c(1.5,10,3))
	layout.show(3)

	par(mar=c(0,0,0,0), ps=8)
	  plot(0,0,xlim=c(0,0.1),ylim=c(0,1),type="n",axes=FALSE,frame.plot=FALSE)
	  text(0.05,0.2,paste0(subset, " Signifcant snps from sigzones of ", trait, sep=""),cex=1.5)
	par(mar=c(0,2,0.7,0.5),ps=8,mgp=c(0.6,0.3,0))
	  plot(res$manhattan,-log10(res$pvalue),xlab="",ylab="",frame.plot=FALSE,xlim=c(deb_emmax,fin_emmax),axes=FALSE, cex =0.3 ,  col=res$col ,pch=16)

	  # label fdr
	  #points(fdr_sub$manhattan,-log10(fdr_sub$pvalue), col="black", cex=1, pch=16)

	  # label sig zones, note that the code for sigzones is a little wonky, so just plot the points in red maybe
	  points(SNP_signif_log4$manhattan,-log10(SNP_signif_log4$pval), col="red", cex=1, pch=16)

	  axis(2,cex.axis=1,lwd=0.5,tck=0.03)
	  plot.window(xlim=c(deb_emmax,fin_emmax),ylim=c(0,max(max(res$log_pval,max(res$log_pval))*1.05)))
	  #lines(res$manhattan,res$lindley,pch=16,col="red",lwd=1)
	  #for(chrom in unique(sigZones$CHROM)){
			#segments(x0=min(res$manhattan[res$chr==chrom]) ,y0= unique(res$th[res$chr==chrom]),
	     #  		 x1=max(res$manhattan[res$chr==chrom]),y1=unique(res$th[res$chr==chrom]),
	      #	     lty=3, col="red")
			#	}  

	box(lwd=0.5)
	axis(4,col="red",lwd=0.5,tck=0.03)
	mtext("-log10(pval)", side=2, line=1, cex=0.8)

	dev.off()



	rm(res)
	rm(list=ls(pattern="cmh")) #test first 
	gc()

}


	
	
