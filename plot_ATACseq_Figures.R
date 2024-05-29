setwd("/genomics/chiragProjects/cellLineProject/atacSeq")

############################################################################################
# Plot to show the distribution of ATAC peaks from each replicate
par (mfrow=c(1,2), mar=c(5,6,5,1) )
barplot( cbind(c(218004,221137,241784, 134320), c(247218,271039,242872,139179) ), beside=T, las=1, names=c("rep1","rep2","rep3","consensus","rep1","rep2","rep3","consensus"),cex.names=0.85, ylab="" )
mtext(side=3, line=1,"HCC1395BL", adj=0.65)
mtext(side=3, line=1,"HCC1395", adj=0.15)
mtext(side=2, line=5,"ATAC peaks")


############################################################################################
# Plot to show distribution of ATAC peaks at genome
A = read.table("outputSummaryGenomicRegionsATACPeaks", header=FALSE)
pdf (file="FigGenomicRegions.pdf", useDingbats=FALSE)
par (mfrow=c(1,2))

plot=barplot( as.matrix(A[,c(2)] ), beside=T, las=1, names=c("HCC1935"), ylim=c(0,80000),col=c("red","blue","grey") )
text(x=plot, y=A$V2, label=round(A$V4,2), pos=3)
legend ("topleft", legend=c("Promoter","Intergenic","Intragenic"), fill=c("red","blue","grey") )

plot=barplot( as.matrix(A[,c(5)] ), beside=T, las=1, names=c("HCC1935BL"), ylim=c(0,80000),col=c("red","blue","grey") )
text(x=plot, y=A$V5, label=round(A$V7,2), pos=3)
legend ("topleft", legend=c("Promoter","Intergenic","Intragenic"), fill=c("red","blue","grey") )
dev.off()


############################################################################################
# Plot to show overlap of ATAC peaks with CGI

A = read.table("outputSummaryCGIATACPeaks", header=FALSE)
pdf (file="FigCGIDistribution.pdf", useDingbats=FALSE)
par (mfrow=c(1,2))
ylab="Frequency"
plot=barplot( as.matrix(A[,c(4)] ), beside=T, las=1, names=c("HCC1935"), ylim=c(0,100), col=c("red","blue"), ylab=ylab  )
text(x=plot, y=A$V4, label=round(A$V4,2), pos=3)
legend ("topleft", legend=c("overlap CGI","nonCGI"), fill=c("red","blue"),cex=1.2 )

plot=barplot( as.matrix(A[,c(7)]
 ), beside=T, las=1, names=c("HCC1935BL"), ylim=c(0,100),col=c("red","blue"), ylab=ylab )
text(x=plot, y=A$V7, label=round(A$V7,2), pos=3)
legend ("topleft", legend=c("overlap CGI","nonCGI"), fill=c("red","blue"),cex=1.2 )

dev.off()



############################################################################################
# Plot to show ATAC peak signal intensity at CGI and nonCGI promtoers

setwd("/genomics/chiragProjects/cellLineProject/atacSeq")

A = read.table("matrixExpHCC1395", header=FALSE)
A$means=rowMeans(A[4:6])
A1_PromoterCGI=  A[A$V2 == "Promoter" & A$V3 == "OverlapCGI",]
A1_PromoterNoCGI=A[A$V2 == "Promoter" & A$V3 == "NonCGI",]
A1_IntragenicCGI   =A[A$V2 == "Intragenic" & A$V3 == "OverlapCGI",]
A1_IntragenicNoCGI =A[A$V2 == "Intragenic" & A$V3 == "NonCGI",]
A1_IntergenicCGI   =A[A$V2 == "Intergenic" & A$V3 == "OverlapCGI",]
A1_IntergenicNoCGI =A[A$V2 == "Intergenic" & A$V3 == "NonCGI",]


pdf (file="FigATACPeakExpLevel.pdf", useDingbats=FALSE, width=14 )
par (mfrow=c(1,2), mar=c(5,6,3,1) )

ylim=c(0,15)
ylab="Expression level- RPKM"

# Define gradinets
library(RColorBrewer)
brewer.pal(n = 8, name = "Reds")
brewer.pal(n = 8, name = "Blues")
col2=c("#084594",  "#6BAED6", "#DEEBF7")
col1=c("#99000D", "#FB6A4A", "#FEE0D2" )


names=c("CGI","CGI","CGI","nonCGI","nonCGI","nonCGI")
boxplot (log2(A1_PromoterCGI$means+1), log2(A1_IntragenicCGI$means+1), log2(A1_IntergenicCGI$means+1), log2(A1_PromoterNoCGI$means+1), log2(A1_IntragenicNoCGI$means+1), log2(A1_IntergenicNoCGI$means+1), ylim=ylim, col=col1, las=1, cex=0.5, notch=T, at=c(1,2,3,5,6,7), names=names,main="HCC1395" )
legend("topright", legend=c("Promoter","Intragenic","Intergenic"), fill=col1)

s1=t.test(log2(A1_PromoterCGI$means+1), log2(A1_PromoterNoCGI$means+1) )
s2=t.test ( log2(A1_IntragenicCGI$means+1), log2(A1_IntragenicNoCGI$means+1) )
s3=t.test ( log2(A1_IntergenicCGI$means+1), log2(A1_IntergenicNoCGI$means+1) )


#boxplot (log2(A1_PromoterCGI$means+1), log2(A1_PromoterNoCGI$means+1) , ylim=ylim, las=1, main="Promoter", names=c("CGI","noncGI"), col=col, ylab=ylab )
#s1=t.test( log2(A1_PromoterCGI$means+1), log2(A1_PromoterNoCGI$means+1) )
#mtext( s1$p.value, side=3, line=-5 )
#boxplot (log2(A1_IntragenicCGI$means+1), log2(A1_IntragenicNoCGI$means+1), ylim=ylim, las=1, main="Intragenic", names=c("CGI","noncGI"), col=col)
#s1=t.test( log2(A1_IntragenicCGI$means+1), log2(A1_IntragenicNoCGI$means+1) )
#mtext( s1$p.value, side=3, line=-5 )
#boxplot (log2(A1_IntergenicCGI$means+1), log2(A1_IntergenicNoCGI$means+1), ylim=ylim, las=1, main="Intergenic", names=c("CGI","noncGI"), col=col)
#1=t.test( log2(A1_IntergenicCGI$means+1), log2(A1_IntergenicNoCGI$means+1) )
#mtext( s1$p.value, side=3, line=-5 )
#dev.off()

A = read.table("matrixExpHCC1395BL", header=FALSE)
A$means=rowMeans(A[4:6])
A1_PromoterCGI=  A[A$V2 == "Promoter" & A$V3 == "OverlapCGI",]
A1_PromoterNoCGI=A[A$V2 == "Promoter" & A$V3 == "NonCGI",]
A1_IntragenicCGI   =A[A$V2 == "Intragenic" & A$V3 == "OverlapCGI",]
A1_IntragenicNoCGI =A[A$V2 == "Intragenic" & A$V3 == "NonCGI",]
A1_IntergenicCGI   =A[A$V2 == "Intergenic" & A$V3 == "OverlapCGI",]
A1_IntergenicNoCGI =A[A$V2 == "Intergenic" & A$V3 == "NonCGI",]

boxplot (log2(A1_PromoterCGI$means+1), log2(A1_IntragenicCGI$means+1), log2(A1_IntergenicCGI$means+1), log2(A1_PromoterNoCGI$means+1), log2(A1_IntragenicNoCGI$means+1), log2(A1_IntergenicNoCGI$means+1), ylim=ylim, col=col2, las=1, cex=0.5, notch=T, at=c(1,2,3,5,6,7), names=names, main="HCC1395BL" )
legend("topright", legend=c("Promoter","Intragenic","Intergenic"), fill=col1)

s1=t.test(log2(A1_PromoterCGI$means+1), log2(A1_PromoterNoCGI$means+1) )
s2=t.test ( log2(A1_IntragenicCGI$means+1), log2(A1_IntragenicNoCGI$means+1) )
s3=t.test ( log2(A1_IntergenicCGI$means+1), log2(A1_IntergenicNoCGI$means+1) )

dev.off()



