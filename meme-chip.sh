#!/bin/bash

# Load MEME-suite module
module load MEME-suite/5.4.1
module load bedtools

# Input a list of bed file names and a list of output file names
bed_file=($1)
outputFile_names=($2)

# Create a new directory for the resulting fasta files
mkdir -p "bed2fasta_files"

# Map the chromosome coordinates in the bed file to the corresponding nucleotide regions in the mouse genome
for ((i=0; i<"${#bed_file[@]}"; i++)); do 
    echo "Running: bed2fasta -o ${bed_file[$i]} mm10.fa ${outputFile_names[$i]}"
    touch "${outputFile_names[$i]}"
    bedtools getfasta -fi "mm10.fa" -bed ${bed_file[$i]} -s -fo ${outputFile_names[$i]}
    mv ${outputFile_names[$i]} "bed2fasta_files"

done

# Run MEME-chip on each of the fasta sequences
cd "bed2fasta_files"
for fasta_file in "${outputFile_names[@]}"; do
    meme-chip $fasta_file

done