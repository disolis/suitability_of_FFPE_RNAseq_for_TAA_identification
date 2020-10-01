# suitability_of_FFPE_RNAseq_for_TAA_identification

## The repository contains some R-code used when writing a master's thesis.


  **Correlation_and_DE_analyses**: 
    correlation (intra- and inter-group variability) and 
    differential expression analysis (FF_od versus FF_control AND FFPE_od versus FFPE_control)


  **LibDiversity**:
    Library diversity between FFPE/FF matching samples:
      Fisher test applying for mapping statistics, duplication rates and read distribution;
      Building plot for inner distance of matching samples;
      Detection of number of all detected protein-coding (PC) genes and number of PC genes consuming 25% of sequencing effort per sample


  **SNVs_calling**:
  * Global (applied for whole transcriptom)
      Finding allele genotypes across all the samples and finding the number of 12 possible substitutions 
      (C>T, G>A, C>A, G>T, C>G, G>C, A>C, T>G, A>G, T>C, A>T, T>A) detected in each sample with respect to the reference genome
  * Local (applied only to the regions of potential tumor-associated antigens)
      Finding the number of all artefacts for the region of TAAs (Note: The study focuses on C> T and G> A)
    
  **STARvsHISAT2**:
  Comparison of two different alignment tools STAR and HISAT2 based on their performance (alignment rate)
 
    
 
