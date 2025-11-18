sig_sl_chr_8=function(lind,POS, th){
	mist = lind
	auxpos = POS
	th = unique(th)
	zones=c(0,0,0)	
	M_loc=which.max(mist)
	if(length(which(mist[1:M_loc]==0))==0){ #the peak is at the beginning of the chrom
			m_loc=1
			zones=rbind(zones, c(auxpos[m_loc],auxpos[M_loc],max(mist)))
			tmp=which.min(which(mist[M_loc+1:length(mist)]==0))
			mist=mist[tmp:length(mist)]
			auxpos=ch1$POS[tmp:length(mist)]
			zones=matrix(zones, ncol=3)
			zones=data.table(beg=zones[,1],end=zones[,2],peak=zones[,3])
		if (nrow(zones)>1){zones=zones[-c(1),]}
	}		

	## split from first peak to run
	x <- which(mist==0)[1]
	ch3 <- ch1[x:nrow(ch1) ,]

	mist3 = ch3$lindley
	auxpos3 = ch3$POS
	zones3=c(0,0,0)	
	M_loc3=which.max(mist3)

	while(max(mist3)>=th){
	  M_loc3=which.max(mist3)
		if(length(which(mist3[1:M_loc3]==0))==0){ #the peak is at the beginning of the chrom
			m_loc3=1
			zones3=rbind(zones3, c(auxpos3[m_loc3],auxpos3[M_loc3],max(mist3)))
			tmp3=which.min(which(mist3[M_loc3+1:length(mist3)]==0))
	  
			mist3=mist3[tmp3:length(mist3)]
			auxpos3=ch3$POS[tmp3:length(mist3)]
			}else{	
				m_loc3=max(which(mist3[1:M_loc3]==0))			
				max3=max(mist3)
				zones3=rbind(zones3, c(auxpos3[m_loc3+1],auxpos3[M_loc3],max3))
				tmp3=which(mist3[M_loc3:length(mist3)]==0) #first 0 score after peak
				if (length(tmp3)>0){
				  auxpos3=auxpos3[c(1:m_loc3,(min(tmp3)+M_loc3):length(mist3))]
				  mist3=mist3[c(1:m_loc3, (min(tmp3)+M_loc3):length(mist3))]
				  }else{ #if the peak is at the end of the chromosome
				    auxpos3=auxpos3[1:m_loc3]
				    mist3=mist3[1:m_loc3]
				    }				
				}
	  }

	zones3=matrix(zones3, ncol=3)
	zones3=data.table(beg=zones3[,1],end=zones3[,2],peak=zones3[,3])
	if (nrow(zones3)>1){zones3=zones3[-1,]}
	if (class(zones)[1]=="numeric") {zlist <- list(zones3)
		}else{zlist <- list(zones, zones3)}

	zones <- rbindlist(zlist)   
    return(zones)
}