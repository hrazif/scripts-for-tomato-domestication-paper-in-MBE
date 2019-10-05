dadi_file=read.table("/home/caicedo/Hamid/20190311_redoing_analyses_for_the_paper/7_dadi_dating/1_full_genome_SP_NECU_and_SLCSLL/same_sample_sizes/no_gene_flow/full_genome.dadi", header=T, check.names=F)


chr_nums=c("01","02","03","04","05","06","07","08","09","10","11","12")

#chunk number
h=0
for (i in 1:length(chr_nums)) {
	temp=dadi_file[dadi_file$Gene==paste0("SL2.50ch",chr_nums[i]),]
	chr_end=temp[nrow(temp),]$Postion
	blocks=seq(1,chr_end, 10000000)
		
	for (j in 1:length(blocks)) {	
		temp2=temp[temp$Postion > blocks[j] & temp$Postion < (blocks[j]+10000000), ]
		h=h+1
		write.table(temp2, paste0("chunk",h,".dadi"), col.names=T, row.names=F, quote=F, sep="\t")
		
	}	
}	