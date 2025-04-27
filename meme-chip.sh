#!/bin/bash

# Load MEME-suite, bedtools and samtools
module load MEME-suite/5.4.1
module load bedtools
module load samtools

# Input a list of bed file names and a list of output file names (do not include .bed or .fa extensions)
bed_file=($1)
outputFile_names=($2)

# Create a .fai file for the mouse genome (mm10.fa)
samtools faidx "mm10.fa"

# Run MEME-chip
for ((i=0; i<"${#bed_file[@]}"; i++)); do 

    # Create output file and directory
    touch "${outputFile_names[$i]}.fa"
    new_dir="meme_${outputFile_names[$i]}"
    mkdir -p $new_dir

    # Map the chromosome coordinates in the bed file to the corresponding nucleotide regions in the mouse genome
    bedtools getfasta -fi "mm10.fa" -bed "${bed_file[$i]}.bed" -s -fo "${outputFile_names[$i]}.fa"
    mv "${outputFile_names[$i]}.fa" $new_dir

    # Run MEME-chip on the generated fasta file
    cd $new_dir
    meme-chip "${outputFile_names[$i]}.fa"
    cd ..

done