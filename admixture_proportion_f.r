

freq_table=data.frame()	
for (i in c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")){

	temp = read.table(paste0("chr",i,".derFreq.tsv.gz"), header=T, as.is=T)
	freq_table=rbind(freq_table, temp)
}



#function for f
f.stat <- function(p1, p2, p3a, p3b) {
    ABBA_numerator <- na.omit((1 - p1) * p2 * p3a)
    BABA_numerator <- na.omit(p1 * (1 - p2) * p3a)

    ABBA_denominator <- na.omit((1 - p1) * p3b * p3a)
    BABA_denominator <- na.omit(p1 * (1 - p3b) * p3a)

    (sum(ABBA_numerator) - sum(BABA_numerator)) /
    (sum(ABBA_denominator) - sum(BABA_denominator))
    }



#populations
P1 <- "SLC_MEX_CA_NSA"
P2 <- "SLC_PER"
P3a <- "SP_PER1"
P3b <- "SP_PER2"



f <- f.stat(freq_table[,P1], freq_table[,P2], freq_table[,P3a], freq_table[,P3b])


#block jackknife
source("./genomics_general/jackknife.R")

#to partition the dataset into blocks on 5MB
block_indices <- get_block_indices(block_size=5e6,
                                   positions=freq_table$position,
                                   chromosomes=freq_table$scaffold)
n_blocks <- length(block_indices)

f_sd <- get_jackknife_sd(block_indices=block_indices,
                         FUN=f.stat,
                         freq_table[,P1], freq_table[,P2], freq_table[,P3a], freq_table[,P3b])


						 
						 
f_err <- f_sd/sqrt(n_blocks)


f_Z <- f / f_err

#95% confidence interval
f_CI_lower <- f - 1.96*f_err
f_CI_upper <- f + 1.96*f_err


c=data.frame("genome-wide f"=f, "Jackknife standard error of f"=f_err, "Z for f"=f_Z)
write.table(c, "stats_output.txt", sep="\t", quote=F, row.names = F)

							   