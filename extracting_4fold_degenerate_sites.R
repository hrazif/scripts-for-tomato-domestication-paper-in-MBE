#@!source("https://bioconductor.org/biocLite.R")
#@!biocLite()
#@!biocLite(c(
#@!  "Biostrings",
#@! "GenomicRanges",
#@! "BSgenome", 
#@!  "rtracklayer",
#@!  "motifRG"
#@!))

#@!install.packages("dplyr")


#@!vignette(package="Biostrings")
#@!vignette("BiostringsQuickOverview", package="Biostrings")


library(Biostrings)     # Provides DNAString, DNAStringSet, etc
library(BSgenome)       # Provides getSeq()
library(GenomicRanges)  # Provides GRanges, etc
library(rtracklayer)    # Provides import() and export()

library(dplyr)          # Pipe %>%



#to load the features file
features <- import.gff3("/home/FrankLab/genomes/tomato/ITAG4.0_gene_models.gff")


# Optional: just retain the columns of metadata we need
mcols(features) <- mcols(features)[,c("type","ID", "Parent")]



#chromosomes
chs=c("ch01", "ch02", "ch03", "ch04", "ch05", "ch06", "ch07", "ch08", "ch09", "ch10", "ch11", "ch12")

df_all=vector()
for(h in 1:12) {
	codon_position=vector()
	#to load each chromosome and convert into a string
	file <- readDNAStringSet(paste("SL4.0", chs[h], ".fa", sep = ""))
	dnastring=paste(file[[1]], collapse = "")

	#to extract CDSs on the plus strand of each chromosome
	cds.chr.plus <- subset(features, type == "CDS" & seqnames==paste("SL4.0", chs[h], sep="") & strand=="+")
	
	mRNA=unique(as.vector(cds.chr.plus$Parent@unlistData))
	for (i in 1: length(mRNA)) {
		mergedCDS=subset(cds.chr.plus, Parent == paste(mRNA[i]))
		#to order CDSs by "start"
		mergedCDS=mergedCDS[with(mergedCDS, order(start))]
		CDS.sequences = Views(dnastring, start=mergedCDS@ranges@start, end=mergedCDS@ranges@start+mergedCDS@ranges@width-1)

		CDS.full.sequence=vector()
		for (j in 1: length(CDS.sequences)) {
			temp.seq=paste(CDS.sequences[[j]])
			CDS.full.sequence=c(CDS.full.sequence,temp.seq)
		}  

		codon.plus=c("GCT","GCC","GCA","GCG","CGT","CGC","CGA","CGG","GGT","GGC","GGA","GGG","CTT","CTC","CTA","CTG","CCT","CCC","CCA","CCG","TCT","TCC","TCA","TCG","ACT","ACC","ACA","ACG","GTT","GTC","GTA","GTG")

		for(k in 1:length(codon.plus)) {
			result <- matchPattern(codon.plus[k], paste(CDS.full.sequence, collapse = ""))
			result=result[(result@ranges@start +2) %% 3 == 0]
			temp=vector()
			result_CDS_1=result[result@ranges@start +2 <= CDS.sequences@ranges@width[1],] 
			temp_1=(CDS.sequences@ranges@start[1]-1 + result_CDS_1@ranges@start+ 2) 
			codon_position=c(codon_position, temp_1)
  			for (l in 2: length(CDS.sequences)) {
				result_CDS_rest=result[result@ranges@start +2 > cumsum(CDS.sequences@ranges@width)[l-1] & result@ranges@start + 2 <= cumsum(CDS.sequences@ranges@width)[l],] 
				temp=(CDS.sequences@ranges@start[l]-1 + result_CDS_rest@ranges@start - cumsum(CDS.sequences@ranges@width)[l-1]) + 2
				codon_position=c(codon_position, temp)
			}
  
		}

	}


	#to extract CDSs on the minus strand of each chromosome
	cds.chr.minus <- subset(features, type == "CDS" & seqnames==paste("SL4.0", chs[h], sep="") & strand=="-")

	mRNA=unique(as.vector(cds.chr.minus$Parent@unlistData))

	for (i in 1: length(mRNA)) {
  
		mergedCDS=subset(cds.chr.minus, Parent == paste(mRNA[i]))
		#to order CDSs by "start"
		mergedCDS=mergedCDS[with(mergedCDS, order(start))]
		CDS.sequences = Views(dnastring, start=mergedCDS@ranges@start, end=mergedCDS@ranges@start+mergedCDS@ranges@width-1)
    
		CDS.full.sequence=vector()
		for (j in 1: length(CDS.sequences)) {
			temp.seq=paste(CDS.sequences[[j]])
			CDS.full.sequence=c(CDS.full.sequence,temp.seq)
		}  
  
		codon.minus=c("CAC","TAC","GAC","AAC","CGT","TGT","GGT","AGT","CGA","TGA","GGA","AGA","CGG","TGG","GGG","AGG","CAG","TAG","GAG","AAG","CCC","TCC","GCC","ACC","CCG","TCG","GCG","ACG","CGC","TGC","GGC","AGC")
  
		for(k in 1:length(codon.minus)) {
			result <- matchPattern(codon.minus[k], paste(CDS.full.sequence, collapse = ""))
			result=result[(result@ranges@start +2) %% 3 == 0]
			temp=vector()
			result_CDS_1=result[result@ranges@start <= CDS.sequences@ranges@width[1],] 
			temp_1=(CDS.sequences@ranges@start[1]-1 + result_CDS_1@ranges@start) 
			codon_position=c(codon_position, temp_1)
    		for (l in 2: length(CDS.sequences)) {
				result_CDS_rest=result[result@ranges@start > cumsum(CDS.sequences@ranges@width)[l-1] & result@ranges@start <= cumsum(CDS.sequences@ranges@width)[l],] 
				temp=(CDS.sequences@ranges@start[l]-1 + result_CDS_rest@ranges@start - cumsum(CDS.sequences@ranges@width)[l-1])
				codon_position=c(codon_position, temp)
			}
    
		}
  
	}

df_chr=data.frame(chrom=paste("SL4.0", chs[h], sep = ""), codon_pos=codon_position)
df_chr= df_chr[with(df_chr, order(codon_pos)), ]
df_all=rbind(df_all, df_chr)


}

write.table(df_all, "chrs_codon_position.txt", quote = F, row.names = F, col.names = F, sep="\t")

