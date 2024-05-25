#! /bin/bash

# Kallisto was used to quantify the expression of transcripts

kallistoIndex=/genomics/Ref_genome/Human/hg38/hg38.ensGene.kallisto.index

HCC1395_rep1_R1=/genomics/1_Projects/FDA_QC/Fastq/bulk_RNA/merge/HCC1395_1/HCC1395_1_R1.fastq.gz
HCC1395_rep1_R2=/genomics/1_Projects/FDA_QC/Fastq/bulk_RNA/merge/HCC1395_1/HCC1395_1_R2.fastq.gz
HCC1395_rep2_R1=/genomics/1_Projects/FDA_QC/Fastq/bulk_RNA/merge/HCC1395_2/HCC1395_2_R1.fastq.gz
HCC1395_rep2_R2=/genomics/1_Projects/FDA_QC/Fastq/bulk_RNA/merge/HCC1395_2/HCC1395_2_R2.fastq.gz
HCC1395_rep3_R1=/genomics/1_Projects/FDA_QC/Fastq/bulk_RNA/merge/HCC1395_3/HCC1395_3_R1.fastq.gz
HCC1395_rep3_R2=/genomics/1_Projects/FDA_QC/Fastq/bulk_RNA/merge/HCC1395_3/HCC1395_3_R2.fastq.gz

HCC1395BL_rep1_R1=/genomics/1_Projects/FDA_QC/Fastq/bulk_RNA/merge/HCC1395BL_1/HCC1395BL_1_R1.fastq.gz
HCC1395BL_rep1_R2=/genomics/1_Projects/FDA_QC/Fastq/bulk_RNA/merge/HCC1395BL_1/HCC1395BL_1_R2.fastq.gz
HCC1395BL_rep2_R1=/genomics/1_Projects/FDA_QC/Fastq/bulk_RNA/merge/HCC1395BL_2/HCC1395BL_2_R1.fastq.gz
HCC1395BL_rep2_R2=/genomics/1_Projects/FDA_QC/Fastq/bulk_RNA/merge/HCC1395BL_2/HCC1395BL_2_R2.fastq.gz
HCC1395BL_rep3_R1=/genomics/1_Projects/FDA_QC/Fastq/bulk_RNA/merge/HCC1395BL_3/HCC1395BL_3_R1.fastq.gz
HCC1395BL_rep3_R2=/genomics/1_Projects/FDA_QC/Fastq/bulk_RNA/merge/HCC1395BL_3/HCC1395BL_3_R2.fastq.gz


# Kallisto quantification
kallisto quant --index=$kallistoIndex --output-dir=/genomics/chiragProjects/cellLineProject/rnaseq/HCC1395_rep1   --threads=4 --plaintext -bias $HCC1395_rep1_R1 $HCC1395_rep1_R2
kallisto quant --index=$kallistoIndex --output-dir=/genomics/chiragProjects/cellLineProject/rnaseq/HCC1395_rep2   --threads=4 --plaintext -bias $HCC1395_rep2_R1 $HCC1395_rep2_R2
kallisto quant --index=$kallistoIndex --output-dir=/genomics/chiragProjects/cellLineProject/rnaseq/HCC1395_rep3   --threads=4 --plaintext -bias $HCC1395_rep3_R1 $HCC1395_rep3_R2
kallisto quant --index=$kallistoIndex --output-dir=/genomics/chiragProjects/cellLineProject/rnaseq/HCC1395BL_rep1 --threads=4 --plaintext -bias $HCC1395BL_rep1_R1 $HCC1395BL_rep1_R2
kallisto quant --index=$kallistoIndex --output-dir=/genomics/chiragProjects/cellLineProject/rnaseq/HCC1395BL_rep2 --threads=4 --plaintext -bias $HCC1395BL_rep2_R1 $HCC1395BL_rep2_R2
kallisto quant --index=$kallistoIndex --output-dir=/genomics/chiragProjects/cellLineProject/rnaseq/HCC1395BL_rep3 --threads=4 --plaintext -bias $HCC1395BL_rep3_R1 $HCC1395BL_rep3_R2


  
# Make tabular matrix of TPM
# Determine number of transcript expressed at different threshold

paste HCC1395_rep1/abundance.tsv   HCC1395_rep2/abundance.tsv   HCC1395_rep3/abundance.tsv   | awk 'BEGIN {OFS="\t"} { out=$1; for (i=5; i<=NF; i+=5) out=out"\t"$i; print out }' | awk '{ if ( ($2 >=0.5 && $3 >=0.5 && $4 >= 0.5) && ($2 >1 || $3 >1 || $4 > 1) ) { print $0 }}' > expressedTranscript_HCC1395
paste HCC1395BL_rep1/abundance.tsv HCC1395BL_rep2/abundance.tsv HCC1395BL_rep3/abundance.tsv | awk 'BEGIN {OFS="\t"} { out=$1; for (i=5; i<=NF; i+=5) out=out"\t"$i; print out }' | awk '{ if ( ($2 >=0.5 && $3 >=0.5 && $4 >= 0.5) && ($2 >1 || $3 >1 || $4 > 1) ) { print $0 }}' > expressedTranscript_HCC1395BL






