library(som)
library(beanplot)
library(RColorBrewer)
library(pheatmap)
library(dendextend)

setwd("/genomics/chiragProjects/cellLineProject/rnaseq")

A = read.table("expressedGenes_HCC1395Representative", header=FALSE)
A1= A[c(7:9)]
M1 = as.matrix(log2(A1+1))
rownames(M1)=A$V1
# Scale data in the range of 0-1
S1=apply(M1, MARGIN = 1, FUN = function(X) (X - min(X))/diff(range(X)))
S2=t(S1)
S2=na.omit(S2)
pdf (file="FigHeatmap_HCC1395.pdf")
pheatmap(S2,cluster_cols=1, fontsize=8, cutree_rows = 2, show_rownames = FALSE,  show_colnames = TRUE )
dev.off()


A = read.table("expressedGenes_HCC1395BLRepresentative", header=FALSE)
A1= A[c(7:9)]
M1 = as.matrix(log2(A1+1))
rownames(M1)=A$V1
# Scale data in the range of 0-1
S1=apply(M1, MARGIN = 1, FUN = function(X) (X - min(X))/diff(range(X)))
S2=t(S1)
S2=na.omit(S2)
pdf (file="FigHeatmap_HCC1395BL.pdf" )
pheatmap(S2,cluster_cols=1, fontsize=8, cutree_rows = 2, show_rownames = FALSE,  show_colnames = TRUE )
dev.off()
                                                                                                                  

A1 = read.table("S1CGI_Top1_BL",   header=FALSE)
A2 = read.table("S1CGI_Top2_BL",   header=FALSE)
A3 = read.table("S1CGI_Top3_BL",   header=FALSE)
A4 = read.table("S1CGI_Top4_BL",   header=FALSE)
A5 = read.table("S1NOCGI_Top1_BL", header=FALSE)
A6 = read.table("S1NOCGI_Top2_BL", header=FALSE)
A7 = read.table("S1NOCGI_Top3_BL", header=FALSE)
A8 = read.table("S1NOCGI_Top4_BL", header=FALSE)

B1 = read.table("S1CGI_Top1",   header=FALSE)
B2 = read.table("S1CGI_Top2",   header=FALSE)
B3 = read.table("S1CGI_Top3",   header=FALSE)
B4 = read.table("S1CGI_Top4",   header=FALSE)
B5 = read.table("S1NOCGI_Top1", header=FALSE)
B6 = read.table("S1NOCGI_Top2", header=FALSE)
B7 = read.table("S1NOCGI_Top3", header=FALSE)
B8 = read.table("S1NOCGI_Top4", header=FALSE)


pdf (file="FigExpQuartile.pdf")

par (mfrow=c(1,2))
col1=c(rep("#99000D",4), rep("#FB6A4A",4))
col2=c(rep("#084594",4), rep("#4292C6",4))
ylim=c(0,14)
boxplot( log2(B1$V10+1), log2(B2$V10+1), log2(B3$V10+1), log2(B4$V10+1),  log2(B5$V10+1), log2(B6$V10+1), log2(B7$V10+1),log2(B8$V10+1),col=col1, las=1 , cex=0.4, at=c(1,2,3,4,6,7,8,9), notch=T, ylim=ylim )
legend("topright", legend=c("CGI","nonCGI"),fill=c("#99000D","#FB6A4A"))

boxplot( log2(A1$V10+1), log2(A2$V10+1), log2(A3$V10+1), log2(A4$V10+1),  log2(A5$V10+1), log2(A6$V10+1), log2(A7$V10+1),log2(A8$V10+1),col=col2, las=1 , cex=0.4, at=c(1,2,3,4,6,7,8,9), notch=T, ylim=ylim  )
legend("topright", legend=c("CGI","nonCGI"),fill=c("#084594","#4292C6"))

dev.off()








