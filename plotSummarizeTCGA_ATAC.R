
A = read.table ("outputSummarizeTCGA_ATAC_HCC1395", header=FALSE)
A$R1=A$V2*100/A$V3
A$R2=A$V4*100/A$V5
A$R3=A$V6*100/A$V7

A$R4=A$V8*100/A$V9
A$R5=A$V10*100/A$V11
A$R6=A$V12*100/A$V13

A1=A[c(14:16)]
rownames(A1)=A$V1
A1=A1[order(A1$R1, decreasing=TRUE),]

A2=A[c(17:19)]
rownames(A2)=A$V1
A2=A2[order(A2$R4, decreasing=TRUE),]

pdf (file="FigATAC_HCC1395.pdf", useDingbats=FALSE, width=16)

par (mfrow=c(1,2))
barplot(t(as.matrix(A1[1:3]) ), beside=T, las=2, col=c("red","blue","grey"), ylim=c(0,100))
box()
barplot(t(as.matrix(A2[1:3]) ), beside=T, las=2, col=c("red","blue","grey"), ylim=c(0,100))
box()

dev.off()



A = read.table ("outputSummarizeTCGA_ATAC_HCC1395BL", header=FALSE)
A$R1=A$V2*100/A$V3
A$R2=A$V4*100/A$V5
A$R3=A$V6*100/A$V7

A$R4=A$V8*100/A$V9
A$R5=A$V10*100/A$V11
A$R6=A$V12*100/A$V13

A1=A[c(14:16)]
rownames(A1)=A$V1
A1=A1[order(A1$R1, decreasing=TRUE),]

A2=A[c(17:19)]
rownames(A2)=A$V1
A2=A2[order(A2$R4, decreasing=TRUE),]

pdf (file="FigATAC_HCC1395BL.pdf", useDingbats=FALSE, width=16 )
par (mfrow=c(1,2))
barplot(t(as.matrix(A1[1:3]) ), beside=T, las=2, col=c("red","blue","grey"), ylim=c(0,100))
box()
barplot(t(as.matrix(A2[1:3]) ), beside=T, las=2, col=c("red","blue","grey"), ylim=c(0,100))
box()

dev.off()







