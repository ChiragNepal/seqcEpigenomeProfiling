#! /bin/bash

# the script generate a coverage plot using deepTools


# Define data path
geneHCC1395=/genomics/chiragProjects/cellLineProject/rnaseq/expressedGenes_HCC1395Representative 
geneHCC1395BL=/genomics/chiragProjects/cellLineProject/rnaseq/expressedGenes_HCC1395BLRepresentative 
CGI=/genomics/Ref_genome/Human/hg38/CGI

# exclude small RNAs
exclude="chrM\|RNA5-8SN2\|SNOR\|ENST00000427161\|MIR"

# ATAC coverage tracks
covHCC1395=/genomics/chiragProjects/cellLineProject/atacSeq/coverage/mergedHCC1395coverage.bw
covHCC1395BL=/genomics/chiragProjects/cellLineProject/atacSeq/coverage/mergedHCC1395BLcoverage.bw
#col2=c("#084594",  "#6BAED6", "#DEEBF7")
#col1=c("#99000D", "#FB6A4A", "#FEE0D2" )




cat $geneHCC1395 | grep -v -e $exclude | awk 'BEGIN {OFS="\t"} { if ($6 == "+") { print $1,$2-300,$2+300 } else { print $1,$3-300,$3+300 }}' | intersectBed -wa -u -a $CGI -b stdin > cgiPromoterHCC1395
cat $geneHCC1395 | grep -v -e $exclude | awk 'BEGIN {OFS="\t"} { if ($6 == "+") { print $1,$2-300,$2+300 } else { print $1,$3-300,$3+300 }}' | intersectBed -v     -a $CGI -b stdin | intersectBed -wa -u -a stdin -b $geneHCC1395 > cgiIntragenicHCC1395

cat $geneHCC1395BL | grep -v -e $exclude | awk 'BEGIN {OFS="\t"} { if ($6 == "+") { print $1,$2-300,$2+300 } else { print $1,$3-300,$3+300 }}' | intersectBed -wa -u -a $CGI -b stdin > cgiPromoterHCC1395BL
cat $geneHCC1395BL | grep -v -e $exclude | awk 'BEGIN {OFS="\t"} { if ($6 == "+") { print $1,$2-300,$2+300 } else { print $1,$3-300,$3+300 }}' | intersectBed -v     -a $CGI -b stdin | intersectBed -wa -u -a stdin -b $geneHCC1395BL > cgiIntragenicHCC1395BL


computeMatrix scale-regions -S $covHCC1395 -R cgiPromoterHCC1395 cgiIntragenicHCC1395 --regionBodyLength 2000 -a 2000 -b 2000 --binSize 100 -p 8 -out matrix2
plotProfile -m matrix2 -out FigureCoverageCGIHCC1395   --regionsLabel "pCGI" "iCGI" --colors "#99000D" "#FB6A4A" --plotFileFormat png --yMin 0 --yMax 30 --startLabel start --endLabel end

computeMatrix scale-regions -S $covHCC1395BL -R cgiPromoterHCC1395BL cgiIntragenicHCC1395BL --regionBodyLength 2000 -a 2000 -b 2000 --binSize 100 -p 8 -out matrix2
plotProfile -m matrix2 -out FigureCoverageCGIHCC1395BL   --regionsLabel "pCGI" "iCGI" --colors "#084594" "#6BAED6" --plotFileFormat png --yMin 0 --yMax 40 --startLabel start --endLabel end
plotProfile -m matrix2 -out FigureCoverageCGIHCC1395BL_2 --regionsLabel "" ""          --colors "#084594" "#6BAED6" --plotFileFormat png --yMin 0 --yMax 40 --startLabel start --endLabel end

# Combined two cell lines
computeMatrix scale-regions -S $covHCC1395  $covHCC1395BL  -R cgiPromoterHCC1395  --regionBodyLength 2000 -a 2000 -b 2000 --binSize 100 -p 8 -out matrix2
plotProfile -m matrix2 -out FigureCoverageCGICombinedPromoter --samplesLabel "" ""    --colors "#99000D" "#084594"  --plotFileFormat png --yMin 0 --yMax 40 --startLabel start --endLabel end --perGroup

computeMatrix scale-regions -S $covHCC1395  $covHCC1395BL  -R cgiIntragenicHCC1395  --regionBodyLength 2000 -a 2000 -b 2000 --binSize 100 -p 8 -out matrix2
plotProfile -m matrix2 -out FigureCoverageCGICombinedIntragenic --samplesLabel  "" ""         --colors "#99000D" "#084594"  --plotFileFormat png --yMin 0 --yMax 40 --startLabel start --endLabel end --perGroup


# Gene body coverage plot
computeMatrix scale-regions -S $covHCC1395 $covHCC1395BL -R genes_HCC1395_CGI genes_HCC1395_NOCGI genes_HCC1395_NOCGI_2 --regionBodyLength 2000 -a 2000 -b 2000 --binSize 100 -p 8 -out matrix2
plotProfile -m matrix2 -out FigureGeneCoverage --regionsLabel "CGI" "nonCGI" --colors "red" "blue" --plotFileFormat png --yMin 0 0 --yMax 50 80


cat $geneHCC1395 |grep -v -e chrM -e RNA5-8SN2 -e SNOR -e ENST00000427161 -e MIR | awk '{ if ($3 - $2 > 500) { print $0 }}' | grep OverlapCGI > genes_HCC1395_CGI
cat $geneHCC1395 |grep -v -e chrM -e RNA5-8SN2 -e SNOR -e ENST00000427161 -e MIR | awk '{ if ($3 - $2 > 1200) { print $0 }}' | grep -v OverlapCGI > genes_HCC1395_NOCGI
# col2=c("#084594",  "#6BAED6", "#DEEBF7")
# col1=c("#99000D", "#FB6A4A", "#FEE0D2" )
computeMatrix scale-regions -S $covHCC1395 -R genes_HCC1395_CGI genes_HCC1395_NOCGI  --regionBodyLength 2000 -a 2000 -b 2000 --binSize 100 -p 8 -out matrix2
plotProfile -m matrix2 -out FigureGeneCoverage --regionsLabel "" "" --colors "#99000D" "#FB6A4A" --yMax 40


cat $geneHCC1395BL |grep -v -e chrM -e RNA5-8SN2 -e SNOR -e ENST00000427161 -e MIR | awk '{ if ($3 - $2 > 500) { print $0 }}'  | grep    OverlapCGI > genes_HCC1395BL_CGI
cat $geneHCC1395BL |grep -v -e chrM -e RNA5-8SN2 -e SNOR -e ENST00000427161 -e MIR | awk '{ if ($3 - $2 > 1200) { print $0 }}' | grep -v OverlapCGI > genes_HCC1395BL_NOCGI
computeMatrix scale-regions -S $covHCC1395BL -R genes_HCC1395BL_CGI genes_HCC1395BL_NOCGI  --regionBodyLength 2000 -a 2000 -b 2000 --binSize 100 -p 8 -out matrix2
plotProfile -m matrix2 -out FigureGeneCoverage --regionsLabel "" "" --colors "#084594" "#6BAED6" --yMax 80





