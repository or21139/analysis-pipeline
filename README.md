# analysis-pipeline
reproduction of analysis:

download the fastq files and use each set of 3 parts as input to FASTAform.sh script which will combine them into 1 fasta file

perform a BLAST search on each of these fasta files and download the results - full descriptions and matching sequences of genomes chosen within blast.fasta files
also download description table data if desired

use gap_filler.sh with the BLAST output file as an input to ensure each genome is continuous (only required for A and C)

determine translation table required for each sequence (all 2 in this case)

combine original fasta files with their corresponding BLAST output in preparation for MSA
  cat sampleA.fasta  sampleA_BLAST_clean.fasta > sampleA_sequences.fasta
  cat sampleB.fasta  sampleB_BLAST.fasta > sampleB_sequences.fasta
  cat sampleC.fasta  sampleC_BLAST_clean.fasta > sampleC_sequences.fasta
  cat sampleD.fasta  sampleD_BLAST.fasta > sampleD_sequences.fasta

use these files as input for sequence_translate.py to translate into amino acid sequence

use the translated sequences as MSA input using MUSCLE

run the output files through msa_cleaner.sh using sample sequences file and the corresponding msa output as input arguments (in that order)

use the resulting file as input for the tree_formation.py, resulting in 4 nexus files containing the trees

upload nexus files to iTOL if you wish to view 
