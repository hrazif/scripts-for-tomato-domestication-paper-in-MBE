#n: number of individuals
n=c(4,        34,      10,     47,      29,    14,     16,     29,      79,      69)

#S: number of segregating sites; i.e. number of SNPs
S=c(1992683,10282318,2548073,6814319,5230980,1121874,2413585,2115931,3536874,1937776)

#genome size (L) =823833640
SNP_density=S/823833640

pops=c("SP_SECU","SP_PER","SP_NECU","SLC_ECU","SLC_PER","SLC_San_Martin","SLC_MEX_CA_NSA","SLC_MEX", "SLL", "SLL_xAdmixed")




theta_watterson=vector()
for(i in 1: length(pops)){
  theta_watterson[i]=SNP_density[i]/sum(1/1:(n[i]-1))
}

dataset=data.frame(pops, n, S, SNP_density, theta_watterson)

write.table(dataset, "20190311_18_theta_waterson_per_pop.txt", row.names = F, col.names = T, sep="\t", quote = F)


#to plot mean Waterson's D per pop

pdf("20190311_18_theta_waterson_per_pop.pdf", width = 6, height = 6)
par(mfrow=c(1,1), mai=c(2, 0.5, 0.3, 0.3))
plot(dataset[,5], xaxt = "n", xlab="", ylab="", main="average Waterson's theta", pch=16, col="red")
axis(1, at=1:10, labels=dataset[,1], cex.axis=1, las=2)
lines(1:10, dataset[,5], col="blue", lwd=2, type="h")
dev.off()
