#!/bin/bash

# converting windows format to unix format by removing the the newline character from windows
sed -i 's/\r//g' "$1"

# remove spaces after pipe character
sed -i 's/| /|/g' "$1"

# because of the next command it is important that the input file is of the shape file.fasta! if it's .fq or something strange stuff will happen
filename=$(basename "$1" .fasta)

# first awk removes al headers that have been seen before
# second awk changes all non nucleotide characters in the sequences

awk '/^>/ {f=!d[$1]; d[$1]=1} f' $1 |  awk '/^>/ {print; next} {gsub(/[^ACGTURYSWKMBDHVN?acgturyswkmbdhvn-]/, "n"); print}' > $filename\_clean.fasta
