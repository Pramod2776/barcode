#!/bin/bash

# Define input files and parameters
REF="./results/test.fasta"
READS="./results/raw/R2.fastq"
OUT_DIR="./output"

# Create output directory if it doesn't exist
mkdir -p $OUT_DIR

# Map reads to reference sequences using BBSplit
./bbmap/bbsplit.sh  in=$READS ref=$REF basename=$OUT_DIR/%_#.fastq outu=$OUT_DIR/unmapped.fastq

# Initialize a file to store read counts
echo -e "Reference\tRead_Count" > $OUT_DIR/read_counts.txt

# Count reads in each output file
for f in $OUT_DIR/*.fastq; do
    REF_NAME=$(basename $f .fastq)
    READ_COUNT=$(./bbmap/reformat.sh  in=$f | grep "Reads" | awk '{print $2}')
    echo -e "$REF_NAME\t$READ_COUNT" >> $OUT_DIR/read_counts.txt
done
