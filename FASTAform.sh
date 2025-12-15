#!/bin/bash

#iterate this for loop for each file part provided for each sample
for i in $@; do
	file_name=$(head -1 $i | grep -o 'sample.') #extract the sample name 
	sequence=$(sed -E '/\+/q; /^@/d' $i | sed -E '/\+/d' | tr -d '\n' | tr -d ' ') #extract the nucleotide sequence from fastq files
		#checking to see whether a fasta file for each sample already exists, if not, adding the sample title
	if test -e "${file_name}.fasta"; then
		echo -e -n ${sequence} >> "${file_name}.fasta" #adding sequence to fasta file
	else
		echo -e -n ">${file_name}\n${sequence}" >> "${file_name}.fasta" #adding title and sequence to fasta file
	fi
done
