#! /bin/bash

atacPeakHCC1395=/genomics/chiragProjects/cellLineProject/atacSeq/consensus_atac_HCC1395
atacPeakHCC1395BL=/genomics/chiragProjects/cellLineProject/atacSeq/consensus_atac_HCC1395BL

output1=outputSummarizeTCGA_ATAC_HCC1395
output2=outputSummarizeTCGA_ATAC_HCC1395BL
rm -rf $output1 $output2

for name in $(ls /genomics/chiragProjects/data/tcga/atac_hg38/* )
do

	ID=$(echo $name | sed 's/\//\t/g' |awk '{ print $NF }' | sed 's/_peakCalls_hg38//g')

	#CGI peaks
	count1=$(cat $atacPeakHCC1395 | grep -w OverlapCGI | grep -w Promoter | wc -l)
	count2=$(cat $atacPeakHCC1395 | grep -w OverlapCGI | grep -w Promoter | intersectBed -wa -u -a stdin -b $name | wc -l)
	count3=$(cat $atacPeakHCC1395 | grep -w OverlapCGI | grep -w Intragenic | wc -l)
        count4=$(cat $atacPeakHCC1395 | grep -w OverlapCGI | grep -w Intragenic | intersectBed -wa -u -a stdin -b $name | wc -l)
	count5=$(cat $atacPeakHCC1395 | grep -w OverlapCGI | grep -w Intergenic | wc -l)
        count6=$(cat $atacPeakHCC1395 | grep -w OverlapCGI | grep -w Intergenic | intersectBed -wa -u -a stdin -b $name | wc -l)

 	# Non CGI peaks
	count11=$(cat $atacPeakHCC1395 | grep -wv OverlapCGI | grep -w Promoter | wc -l)
        count12=$(cat $atacPeakHCC1395 | grep -wv OverlapCGI | grep -w Promoter | intersectBed -wa -u -a stdin -b $name | wc -l)
        count13=$(cat $atacPeakHCC1395 | grep -wv OverlapCGI | grep -w Intragenic | wc -l)
        count14=$(cat $atacPeakHCC1395 | grep -wv OverlapCGI | grep -w Intragenic | intersectBed -wa -u -a stdin -b $name | wc -l)
        count15=$(cat $atacPeakHCC1395 | grep -wv OverlapCGI | grep -w Intergenic | wc -l)
        count16=$(cat $atacPeakHCC1395 | grep -wv OverlapCGI | grep -w Intergenic | intersectBed -wa -u -a stdin -b $name | wc -l)

	echo $ID $count2 $count1 $count4 $count3 $count6 $count5 $count12 $count11 $count14 $count13 $count16 $count15 | awk '{ out=$1; for ( i =2; i<= NF; i++) out=out"\t"$i; print out }' >> $output1



        #CGI peaks
        count1=$(cat $atacPeakHCC1395BL | grep -w OverlapCGI | grep -w Promoter | wc -l)
        count2=$(cat $atacPeakHCC1395BL | grep -w OverlapCGI | grep -w Promoter | intersectBed -wa -u -a stdin -b $name | wc -l)
        count3=$(cat $atacPeakHCC1395BL | grep -w OverlapCGI | grep -w Intragenic | wc -l)
        count4=$(cat $atacPeakHCC1395BL | grep -w OverlapCGI | grep -w Intragenic | intersectBed -wa -u -a stdin -b $name | wc -l)
        count5=$(cat $atacPeakHCC1395BL | grep -w OverlapCGI | grep -w Intergenic | wc -l)
        count6=$(cat $atacPeakHCC1395BL | grep -w OverlapCGI | grep -w Intergenic | intersectBed -wa -u -a stdin -b $name | wc -l)

        # Non CGI peaks
        count11=$(cat $atacPeakHCC1395BL | grep -wv OverlapCGI | grep -w Promoter | wc -l)
        count12=$(cat $atacPeakHCC1395BL | grep -wv OverlapCGI | grep -w Promoter | intersectBed -wa -u -a stdin -b $name | wc -l)
        count13=$(cat $atacPeakHCC1395BL | grep -wv OverlapCGI | grep -w Intragenic | wc -l)
        count14=$(cat $atacPeakHCC1395BL | grep -wv OverlapCGI | grep -w Intragenic | intersectBed -wa -u -a stdin -b $name | wc -l)
        count15=$(cat $atacPeakHCC1395BL | grep -wv OverlapCGI | grep -w Intergenic | wc -l)
        count16=$(cat $atacPeakHCC1395BL | grep -wv OverlapCGI | grep -w Intergenic | intersectBed -wa -u -a stdin -b $name | wc -l)

        echo $ID $count2 $count1 $count4 $count3 $count6 $count5 $count12 $count11 $count14 $count13 $count16 $count15 | awk '{ out=$1; for ( i =2; i<= NF; i++) out=out"\t"$i; print out }' >> $output2

done


