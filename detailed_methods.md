## Detailed methods for standard riboCleaner workflow

### Identify false gene models with rDNA sequences

**a.**	barrnap (version 0.9, https://github.com/tseemann/barrnap) is used to determine rDNA sequences that code rRNAs (28s, 18s and 5.8s) from the reference genome of interest.
  
**b.**	Concatenate nearby (5 kb) rRNA genes to determine rDNA regions in the reference genome (including spacer sequences) using BEDtools.
  
**c.**	Overlap identified rDNA regions with gene models from the reference genome using BEDtoods. A gene model is considered to be a false gene model if it overlaps with an rDNA region.
  
**d.**	Use mRNA sequences of all genes from the reference genome as a query to BLAST against false gene model sequences, to obtain additional false gene models with rDNA contamination. A query mRNA sequence is considered to be a false model if it has a BLAST hit to a sequence from c that covers >= 90% of the query sequence at >= 90% identity.
  
**e.**	Merge all false gene models from **c** and **d** together
  
### Quantify rRNA contamination

**a.**	Use salmon to map each input RNAseq library to mRNA sequences from the reference genome and calculate TPM values for reference gene models (including false gene models). 

**b.**	Calculate total value and percentage of TPM from the collection of false gene models by adding together the values from each false gene model and then dividing by total TPM of the library

**c.**	Adjust TPM values of all other gene models by rescaling to a million total reads per library after removing the false gene models

### Visualization of the rRNA contamination results

**a.**	Generate the overall rRNA contamination percentage pie chart

**b.**	Generate bar plot of rRNA contamination percentage for each RNAseq library




