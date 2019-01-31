#!/bin/bash

#authored by Colin.
#08/10/2018

set -e
echo "PROCESSING"


READ1="MW1_ATCACG_L007_R1_001.fastq.gz,MW2_CGATGT_L007_R1_001.fastq.gz,MW3_TTAGGC_L007_R1_001.fastq.gz,MW4_TGACCA_L007_R1_001.fastq.gz,MW8_CCGTCC_L007_R1_001.fastq.gz,MW7_TAGCTT_L007_R1_001.fastq.gz,MW6_GATCAG_L007_R1_001.fastq.gz,MW5_ACTTGA_L007_R1_001.fastq.gz"
READ2="MW1_ATCACG_L007_R2_001.fastq.gz,MW2_CGATGT_L007_R2_001.fastq.gz,MW3_TTAGGC_L007_R2_001.fastq.gz,MW4_TGACCA_L007_R2_001.fastq.gz,MW8_CCGTCC_L007_R2_001.fastq.gz,MW7_TAGCTT_L007_R2_001.fastq.gz,MW6_GATCAG_L007_R2_001.fastq.gz,MW5_ACTTGA_L007_R2_001.fastq.gz"


echo $READ1
echo $READ2	

while read LINE
do	
		hisat2 -q -x CSFBv1.2 -1 $READ1 -2 $READ2 -S $LINE.sam
		hisat2 -q -x CSFBv1.2 -1 $READ1 -2 $READ2 -S $LINE.sam
		hisat2 -q -x CSFBv1.2 -1 $READ1 -2 $READ2 -S $LINE.sam
		hisat2 -q -x CSFBv1.2 -1 $READ1 -2 $READ2 -S $LINE.sam
		hisat2 -q -x CSFBv1.2 -1 $READ1 -2 $READ2 -S $LINE.sam
		hisat2 -q -x CSFBv1.2 -1 $READ1 -2 $READ2 -S $LINE.sam
		echo "hisat2 -q -x CSFBv1.2 -1 $READ1 -2 $READ2 -S $LINE.sam"
done <$1
#must output individual commands because hisat2 cannot recognize MW1_ATCACG_L007_R1_001_R1 and MW2_CGATGT_L007_R1_001 as two Read1 files
#but will think of them as the same file.

samtools view  -bS $LINE.sam > $LINE.bam
samtools sort $LINE.bam $LINE.sort

featureCounts -a CSFBv1.5.gff -o counts.txt $LINE.sort.bam 

echo "COMPLETE"

exit 


#if
#	[-e $LINE.bam]
#then
#	featureCounts -a CSFBv1.5.gff -o counts.txt $LINE.bam
#else
#	featureCounts -a CSFBv1.5.gff -o counts.txt $LINE.sam
#fi

#done

#if bam files have been created and exists then use bam files for feature count, else use sam files for feature counts.

#if
#	[-e $LINE.bam]
#then	featureCounts -a CSFBv1.5.gff -o counts.txt $LINE.bam
#
#else
#	echo "ERROR!!!!!! Sam > Bam failed miserably."
#fi



