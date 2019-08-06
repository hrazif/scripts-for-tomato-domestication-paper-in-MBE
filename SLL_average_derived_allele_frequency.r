pop_names=c("SLC_ECU", "SLC_PER", "SLC_San_Martin", "SLC_MEX_CA_NSA", "SLC_MEX", "SLL")

pop_names_simple=c("SLC ECU","SLC PER","SLC San Martin","SLC MEX-CA-NSA","SLC MEX","SLL")


input=read.table("input.dadi", header=T, check.names=F)



#keeping only snps where ancestral and derived alleles are different; to be able to distinguish ancestral and derived alleles. 

temp1=input[input$Ref!=input$OUT,]
temp2=input[input$Ref==input$OUT,]

ancestral_1=temp1[,(length(pop_names)+5):(length(pop_names)+5+length(pop_names)-1)]
ancestral_2=temp2[,4:(length(pop_names)+3)]
ancestral=rbind(ancestral_1,ancestral_2)

derived_1=temp1[,4:(length(pop_names)+3)]
derived_2=temp2[,(length(pop_names)+5):(length(pop_names)+5+length(pop_names)-1)]
derived=rbind(derived_1,derived_2)



#frequencies of derived alleles

freqs_derived=data.frame(matrix(0,nrow(derived),ncol(derived)))

for (i in 1:ncol(derived)){
	freqs_derived[,i]=derived[,i]/(derived[,i]+ancestral[,i])
}
colnames(freqs_derived)=colnames(derived)
freqs_derived=freqs_derived[pop_names]



means_derived=vector()
for (i in 1:ncol(freqs_derived)){

	means_derived[i]=mean(na.omit(freqs_derived[,i]))

}


#frequencies of ancestral alleles

freqs_ancestral=data.frame(matrix(0,nrow(ancestral),ncol(ancestral)))

for (i in 1:ncol(ancestral)){
	freqs_ancestral[,i]=ancestral[,i]/(derived[,i]+ancestral[,i])
}
colnames(freqs_ancestral)=colnames(ancestral)
freqs_ancestral=freqs_ancestral[pop_names]



means_ancestral=vector()
for (i in 1:ncol(freqs_ancestral)){

	means_ancestral[i]=mean(na.omit(freqs_ancestral[,i]))

}



library(vioplot)
pdf("20190311_9_3_violin_derived_allele_freq.pdf", width=6, height = 6)

par(las=1,bty="l")  ## my preferred setting
par(mar=c(8,4.1,4.1,2.1))#sets the bottom, left, top and right margins respectively of the plot region in number of lines of text.

## set up empty plot
y=seq(0,1,0.1)
x=y*10
plot(x,y, axes = F,ann = F, type = "n")
vioplot(freqs_derived[,1],freqs_derived[,2],freqs_derived[,3],freqs_derived[,4],freqs_derived[,5],freqs_derived[,6],freqs_derived[,7],freqs_derived[,8],freqs_derived[,9],add=TRUE, col="red", pchMed=NA)
points(1:length(pop_names), means_derived, col="white", pch=19)
at=1:length(pop_names)
axis(side=1, labels=F, at=at, pos=-.05)
text(cex=1, x=at, y=-0.09, pop_names_simple,xpd=T, srt=45, adj = 1, col=c("gray", "turquoise2", "purple", "lightsalmon1", "darkgreen", "magenta", "orange", "blue", "red"))
axis(side=2, labels=F, at=seq(0,1,0.1), pos=-0.5)
text(cex=1, x=-1, y=c(0,0.5,1), c("0.0","0.5","1.0"), xpd=T, srt=90, adj = 0.5)


dev.off()


library(vioplot)
pdf("20190311_9_3_violin_ancestral_allele_freq.pdf", width=6, height = 6)

par(las=1,bty="l")  ## my preferred setting
par(mar=c(8,4.1,4.1,2.1))#sets the bottom, left, top and right margins respectively of the plot region in number of lines of text.

## set up empty plot
y=seq(0,1,0.1)
x=y*10
plot(x,y, axes = F,ann = F, type = "n")
vioplot(freqs_ancestral[,1],freqs_ancestral[,2],freqs_ancestral[,3],freqs_ancestral[,4],freqs_ancestral[,5],freqs_ancestral[,6],freqs_ancestral[,7],freqs_ancestral[,8],freqs_ancestral[,9],add=TRUE, col="red", pchMed=NA)
points(1:length(pop_names), means_ancestral, col="white", pch=19)
at=1:length(pop_names)
axis(side=1, labels=F, at=at, pos=-.05)
text(cex=1, x=at, y=-0.09, pop_names_simple,xpd=T, srt=45, adj = 1, col=c("gray", "turquoise2", "purple", "lightsalmon1", "darkgreen", "magenta", "orange", "blue", "red"))
axis(side=2, labels=F, at=seq(0,1,0.1), pos=-0.5)
text(cex=1, x=-1, y=c(0,0.5,1), c("0.0","0.5","1.0"), xpd=T, srt=90, adj = 0.5)


dev.off()







