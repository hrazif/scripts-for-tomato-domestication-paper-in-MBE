pops=c("SP_SECU", "SP_PER", "SP_NECU", "SLC_ECU", "SLC_PER", "SLC_San_Martin", "SLC_MEX_CA_NSA", "SLC_MEX", "SLL")
dadi_file=read.table("/home/caicedo/Hamid/20190311_redoing_analyses_for_the_paper/32_redoing_analysis_of_private_alleles/2nd_try_sweep_input/sweeps_input_simple.vcf.data", header=T, check.names=F)


#counts of reference alleles
ref_alleles_counts=dadi_file[,4:(length(pops)+3)]
#to sort columns by populations
ref_alleles_counts=ref_alleles_counts[,pops]


#counts of alternate alleles
alt_alleles_counts=dadi_file[,(length(pops)+5):((length(pops)+5)+length(pops)-1)]
#to sort columns by populations
alt_alleles_counts=alt_alleles_counts[,pops]


#to find private reference alleles

summary_pops_DAF=data.frame(pops=pops, DAF=0, num_private=0)
for (i in 1:length(pops)) {
	
	
	#cases where counts for one population is not zero, but zero for all other popoulations
	temp1=alt_alleles_counts[alt_alleles_counts[,pops[i]]!=0 & apply(alt_alleles_counts[,pops[-i]], 1, sum)==0,]
	
	#to add ref and outgroup information
	temp1$Ref=dadi_file[row.names(temp1),]$Ref
	temp1$OUT=dadi_file[row.names(temp1),]$OUT
	#keeping only the sites where ref and outgroup are the same; thus the ALT is the derived allele
	temp1=temp1[temp1$Ref==temp1$OUT,]
			
	temp2=ref_alleles_counts[ref_alleles_counts[,pops[i]]!=0 & apply(ref_alleles_counts[,pops[-i]], 1, sum)==0,]
	temp2$Ref=dadi_file[row.names(temp2),]$Ref
	temp2$OUT=dadi_file[row.names(temp2),]$OUT
	#keeping only the sites where ref and outgroup are different; thus the REF is the derived allele
	temp2=temp2[temp2$Ref!=temp2$OUT,]
		
	#all private alleles in one pop
	private_alleles_pop=dadi_file[sort(c(rownames(temp1),rownames(temp2))), ]
	
	
	#to sort
	#ref allele counts
	temp3=private_alleles_pop[,4:(length(pops)+3)]
	#to sort columns by populations
	temp3=temp3[,pops]
	#alt allele counts
	temp4=private_alleles_pop[,(length(pops)+5):((length(pops)+5)+length(pops)-1)]
	temp4=temp4[,pops]
	temp5=cbind(temp3[i],temp4[i])
	temp5$Ref=private_alleles_pop$Ref
	temp5$OUT=private_alleles_pop$OUT
	temp5$chr=private_alleles_pop$Gene
	temp5$pop=private_alleles_pop$Postion
	temp5$REF=private_alleles_pop$Allele1
	temp5$ALT=private_alleles_pop$Allele2
		
	temp5$x=0
	temp5[temp5$Ref==temp5$OUT,]$x=temp5[temp5$Ref==temp5$OUT,2]
	temp5[temp5$Ref!=temp5$OUT,]$x=temp5[temp5$Ref!=temp5$OUT,1]
		
	#derived allele frequency (DAF)
	temp5$DAF=temp5$x/(temp5[,1]+temp5[,2])		
	
	summary_pops_DAF$num_private[i]=nrow(temp5)
	summary_pops_DAF$DAF[i]=mean(temp5$DAF)
	
	colnames(temp5)[1:2]=c("ref_allele_count", "alt_allele_count")
	#to check progress
	print(i)
	write.table(temp5, paste0("private_novel_alleles_pop_",pops[i],".txt"), quote=F, row.names=F, col.names=T)
	
}	
	
write.table(summary_pops_DAF, "summary_pops_novel_alleles_DAF.txt", quote=F, row.names=F, col.names=T)	
