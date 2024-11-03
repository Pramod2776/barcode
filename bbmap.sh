./bbmap/bbduk.sh in=./bbmap/raw/R2.fastq literal=GCGAGGAGGCCTTGCGTATA k=24 mm=f hdist=0 out=./results_mapping/unmatched.fq outm=./results_mapping/matched.fq rcomp=f
seqtk seq -L 2 ./results_mapping/matched.fq > ./results_mapping/matched_cleanseqs.fastq
echo $(cat ./results_mapping/matched_cleanseqs.fastq |wc -l)/4|bc > ./results_mapping/matched_cleanseqs_readcount.txt
