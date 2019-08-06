df=read.table("coordinates.txt", header = T, sep="\t")

populations = read.table("populations.txt", header=T)

df=merge(df,populations,by="accession", all=T, sort=F)




df=df[!is.na(df$latitude),]
df=df[!is.na(df$population),]

df$color=0
df$symbol=1
pop_names=c("SP_SECU","SP_PER","SP_NECU","SLC_ECU","SLC_PER","SLC_San_Martin","SLC_MEX_CA_NSA","SLC_MEX","SLL")
pop_colors=c("gray", "turquoise2", "darkmagenta", "lightsalmon1", "darkgreen", "hotpink", "orange", "blue", "red4")
for (i in 1:length(pop_names)) {
df[df$population==pop_names[i],]$color=pop_colors[i]
}

write.table(df, "coordinates_populations.txt", sep="\t", quote = F, row.names = F)
#to write csv
write.table(df, "coordinates_populations.csv", sep=",", quote = F, row.names = F)




library(maps)

pdf("20190311_4_map.pdf", width = 5, height = 4)
par(mai=c(0,0,0,0))
map(database = "world", regions = ".", xlim=c(-110.4,-66.6), ylim=c(-24, 30), fill=T, col="white", resolution=0, bg="transparent", mar = c(0.1, 0.1, 0.1, 0.1))

points(df$longitude, df$latitude, col=df$color, pch=df$symbol, cex=0.5)

#get x and y from google maps by clicking on the country names
text(c(-102, -87, -88.5, -90.2, -85, -83.5, -72.14, -77.0, -74.53), c(24, 15, 13.3, 15.5, 13, 9.3, 4.70, -0.4, -8.5), labels = c("MEX", "HND", "SLV", "GTM", "NIC", "CRI", "COL", "ECU", "PER"), cex=0.4)
x=c("SP SECU","SP PER","SP NECU","SLC ECU","SLC PER","SLC San Martin","SLC MEX-CA-NSA","SLC MEX","SLL", "uncertain")
#legend("bottomleft",x, col=pop_colors, pch=1, cex=0.5 )

dev.off()





