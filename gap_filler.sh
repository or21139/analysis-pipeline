#!/bin/bash

#each input file requires consistent accession codes for each sample BLAST hit

for fasta_file in "${@:1}"; do #repeat the entire script for every input file
	mapfile -t organism_code < <(grep -oP '^>[^:]*' $fasta_file) #create an array of the accession codes
	mapfile -t genome_positions < <(grep -oP ':\K[0-9]+\-[0-9]+' $fasta_file) #create an array of genome position markers
	mapfile -t other_header_info < <(grep -oP '[0-9]+\-[0-9]+\K.*' $fasta_file) #create an array of the rest of each header

	line_count=$(wc -l < $fasta_file) #find how long the file is
	index=0 #start a counter for the index of an array

	#extract each line from the file, if it is a sequence line, append it to the item in the array
	#if a header line, continue to next index in array
	for ((lines=1; lines<=line_count; lines++)); do
		current_line=$(sed -n "${lines}p" $fasta_file)
		if [[ "${current_line:0:1}" == ">" ]]; then
			if [[ lines != 1 ]]; then
				((index++))
			fi
		else
			sequences[${index}-1]+="${current_line}"
		fi
	done

	file_name="${fasta_file%.*}_clean.fasta" #create new file name to output result into

	#by obvserving each accession code, compare whether each sequence is from the same genome as the previous
	for i in "${!organism_code[@]}"; do
		 if [[ $i -eq 0 ]]; then #for the first sequence, there is no previous
                        start_seq=$(grep -oP '^\d+' <<< "${genome_positions[i]}" ) #denotes the starting position of the sequence in the genome
                        end_seq=$(grep -oP '(?<=-)\d+' <<< "${genome_positions[i]}") #denotes the ending position of the sequence in the genome
                        sequence_str="${sequences[i]}" #starts a string to output the final concatonated sequence
		else
			prev=$((i-1)) #defines the previous index
			if [[ "${organism_code[i]}" == "${organism_code[prev]}" ]]; then #compares accession code to previous, if the same:

				#defines genome positions of current and previous sequences for comparison
				prev_start=$(grep -oP '^\d+' <<< "${genome_positions[prev]}" )
				prev_end=$(grep -oP '(?<=-)\d+' <<< "${genome_positions[prev]}")
				current_start=$(grep -oP '^\d+' <<< "${genome_positions[i]}" )
				current_end=$(grep -oP '(?<=-)\d+' <<< "${genome_positions[i]}")

				#depending on whether the BLAST	gets a hit on the forward or backward strand,
				#the sequences will either be in ascending or descending order when downloaded

				if [[ ${prev_end} -lt ${current_start} ]]; then #if ascending
					gap_size=$((${current_start} - ${prev_end}-1)) #calculates the gap between the sequences in the genome
					end_seq=${current_end} #corrects the end position of the sequence
					sequence_str+="$(printf '%.0sN' $(seq 1 $gap_size))${sequences[i]}" #adds Ns to fill the gap along with appending the new sequence

				else #if descending
					gap_size=$((${prev_start} - ${current_end}-1)) #calculates gap size between the sequences in the genome
					start_seq=${current_start} #corrects for the new start position in the sequence
					sequence_str="${sequences[i]}$(printf '%.0sN' $(seq 1 $gap_size))${sequence_str}" #prepends sequence with Ns to fill the gap
				fi

				#if the final accession code in the file, print to file anyway
				if [[ $i -eq $((${#organism_code[@]} - 1)) ]]; then
					echo -e "${organism_code[i]}:${start_seq}-${end_seq}${other_header_info[i]}\n${sequence_str}" >> "$file_name"
				fi
			else #if not the same as previous accession, print out the resulting final sequence with gaps filled
				echo -e "${organism_code[prev]}:${start_seq}-${end_seq}${other_header_info[prev]}\n${sequence_str}" >> "$file_name"
				#reset variables for next iteration
				start_seq=$(grep -oP '^\d+' <<< "${genome_positions[i]}" )
	                       	end_seq=$(grep -oP '(?<=-)\d+' <<< "${genome_positions[i]}")
        	       	        sequence_str="${sequences[i]}"
			fi
		fi
	done
done
