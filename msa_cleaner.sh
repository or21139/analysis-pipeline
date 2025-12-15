#/bin/bash

mapfile -t organism_code < <(grep -oP '^>[^ ]*' $1) #create an array of the accession codes
mapfile -t other_header_info < <(grep -oP '[0-9]+\-[0-9]+\K.*' $1) #create an array of the rest of each header
other_header_info=("" "${other_header_info[@]}")

#add the header info to each header line
for index in ${!organism_code[@]}; do
	sed -i -e "s#${organism_code[index]}#${organism_code[index]}${other_header_info[index]}#" $2
done
