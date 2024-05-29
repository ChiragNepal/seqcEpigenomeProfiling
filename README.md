# seqcEpigenomeProfiling

This repository contains scripts and data for summarizing and plotting epigenome profiles HCC1395 and HCC1395BL cell lines.
This repository contains scripts and data for summarizing and plotting HCC1395 and HCC1395BL cell lines for comparing with TCGA ATAC-seq data.

## Files
- `outputSummaryCGIATACPeaks`: Summary output of ATAC peaks overlapping CGI.
- `outputSummaryGenomicRegionsATACPeaks`: Summary output of ATAC peaks overlapping promoter, intragenic and intergenic regions.
- `outputSummarizeTCGA_ATAC_HCC1395`: Summary output for the HCC1395 cell line.
- `outputSummarizeTCGA_ATAC_HCC1395BL`: Summary output for the HCC1395BL cell line.
- `plotSummarizeTCGA_ATAC.R`: R script for plotting the summarized data.
- `runSummarizeTCGA_ATAC.sh`: Shell script for running the summarization pipeline.


## Usage

1. Run the summarization pipeline using `runSummarizeTCGA_ATAC.sh`.
2. Use `plotSummarizeTCGA_ATAC.R` to generate plots from the summary outputs.
