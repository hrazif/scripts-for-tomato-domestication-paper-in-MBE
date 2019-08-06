

freq_table=data.frame()	
for (i in c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")){

	temp = read.table(paste0("/home/caicedo/Hamid/20180813_redoing_analyses_based_on_new_tree/6_ABBA_BABA_using_Simons_new_scripts/chr",i,".derFreq.tsv.gz"), header=T, as.is=T)
	freq_table=rbind(freq_table, temp)
}



D.stat <- function(p1, p2, p3) {
    ABBA <- na.omit((1 - p1) * p2 * p3)
    BABA <- na.omit(p1 * (1 - p2) * p3)
    (sum(ABBA) - sum(BABA)) / (sum(ABBA) + sum(BABA))
    }	
	

P1 <- "SLC_MEX_CA_NSA"
P2 <- "SLC_PER"
P3 <- "SP_PER"


D <- D.stat(freq_table[,P1], freq_table[,P2], freq_table[,P3])	


#for jackknife
source("./test_with_Simons_scripts/genomics_general/jackknife.R")

#to partition the dataset into blocks on 5MB
block_indices <- get_block_indices(block_size=5e6,
                                   positions=freq_table$position,
                                   chromosomes=freq_table$scaffold)

n_blocks <- length(block_indices)

print(paste("Genome divided into", n_blocks, "blocks."))

D_sd <- get_jackknife_sd(block_indices=block_indices,
                         FUN=D.stat,
                         freq_table[,P1], freq_table[,P2], freq_table[,P3])

print(paste("D standard deviation = ", round(D_sd,4)))


D_err <- D_sd/sqrt(n_blocks)
D_Z <- D / D_err

c=data.frame("genome-wide D"=D, "Jackknife standard error of D"=D_err, "Jackknife Z of D"=D_Z)
write.table(c, "stats_output.txt", sep="\t", quote=F, row.names = F)

								   