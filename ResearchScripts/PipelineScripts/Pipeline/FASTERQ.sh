#!/bin/bash
#first argument: Path to list of all accessions
#Second argument: Output folder for fastq files
#third argument: Path of user repository(set by SRA configuration)


path="$1"



prefetch --option-file $path



while read p; do

	        fasterq-dump $p -O $2 && rm $3/sra/$p.sra
done < "$path"

