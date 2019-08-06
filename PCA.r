library(gdsfmt)
library(SNPRelate)

populations=read.table("populations.txt", header=T)

#plink files; InitFilt_exWilds_exF1_xInd25p_xProblematic_Xinvariant_XSNPs10p

#input plink files also used for TreeMix
bed='/home/caicedo/Hamid/20170417/11_fastSTRUCTURE_like_6_removing_more_problematic_accessions [good]/plink.bed'
fam='/home/caicedo/Hamid/20170417/11_fastSTRUCTURE_like_6_removing_more_problematic_accessions [good]/plink.fam'
bim='/home/caicedo/Hamid/20170417/11_fastSTRUCTURE_like_6_removing_more_problematic_accessions [good]/plink.bim'

###snpgdsBED2GDS(bed, fam, bim, "20180806_1_SNPs.gds")

genofile <- snpgdsOpen("/home/caicedo/Hamid/20180806_redoing_analyses_based_on_new_tree/1_PCA/20180806_1_SNPs.gds")



# Take out snp.id
snp.id=read.gdsn(index.gdsn(genofile, "snp.id"))

# Run PCA
pca <- snpgdsPCA(genofile, snp.id=snp.id, autosome.only=F, remove.monosnp=F, num.thread=10)
head(pca$eigenvect)


# variance proportion (%)
(pc.percent <- pca$varprop*100)
(head(round(pc.percent, 2)))


# make a data.frame
tab <- data.frame(accession = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  EV3 = pca$eigenvect[,3],
				  EV4 = pca$eigenvect[,4],
				  EV5 = pca$eigenvect[,5],
				  stringsAsFactors = FALSE)



tab=merge(tab, populations, by="accession")
tab$population=as.character(tab$population)




tab$color=0
tab$symbol=1
pop_names=c("SP_SECU","SP_PER","SP_NECU","SLC_ECU","SLC_PER","SLC_San_Martin","SLC_MEX_CA_NSA","SLC_MEX","SLL")
pop_colors=c("gray", "turquoise2", "darkmagenta", "lightsalmon1", "darkgreen", "hotpink", "orange", "blue", "red4")

for (i in 1:9) {
tab[tab$population==pop_names[i],]$color=pop_colors[i]
}

write.table(tab, "20190311_1_PCA_output_tab.txt", sep="\t", row.names = F, quote = F)

# Draw mulitple plots with PC1, PC2, etc.
pdf("PCA2-1_20190311_1.pdf", width = 6, height = 6)
plot(tab$EV2, tab$EV1, xlab="PC2", col=tab$color, ylab="PC1", pch=tab$symbol, frame.plot = F)
#legend("topleft",pop_names, col=pop_colors, pch=1, cex=0.5)
dev.off()

pdf("PCA3-2_20190311_1.pdf", width = 6, height = 6)
plot(tab$EV3, tab$EV2, xlab="PC3", col=tab$color, ylab="PC2", pch=tab$symbol, frame.plot = F)
#legend("topleft",pop_names, col=pop_colors, pch=1, cex=0.5)
dev.off()

pdf("PCA4-3_20190311_1.pdf", width = 6, height = 6)
plot(tab$EV4, tab$EV3, xlab="PC4", col=tab$color, ylab="PC3", pch=tab$symbol, frame.plot = F)
#legend("topleft",pop_names, col=pop_colors, pch=1, cex=0.5)
dev.off()

pdf("PCA5-4_20190311_1.pdf", width = 6, height = 6)
plot(tab$EV5, tab$EV4, xlab="PC5", col=tab$color, ylab="PC4", pch=tab$symbol, frame.plot = F)
#legend("topleft",pop_names, col=pop_colors, pch=1, cex=0.5)
dev.off()



