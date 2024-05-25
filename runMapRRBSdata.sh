#! /bin/bash

# This is command used to map the RRBS data.


ref=ref=/genomics/Ref_genome/Human/NCBI_GRCh38/

# Trim sequence with trim galore
trim_galore -phred33 --non_directional --cores 8 --rrbs --paired $FASTQDIR/$pair1 $FASTQDIR/$pair2


# Bismark mapping
bismark --phred33 --bam $ref --non_directional --multicore 12 --unmapped -1 $pair1 -2 $pair2


# Extraction methyaltion call 
# --CX_context also give CH methayltion call on coverage
bismark_methylation_extractor $bamFile --multicore 16 --ignore 2 --ignore_3prime 2 --comprehensive --bedGraph --merge_non_CpG --CX_context --buffer_size 10G


