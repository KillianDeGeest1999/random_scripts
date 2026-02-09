#!/bin/bash

# converting windows format to unix format by removing the the newline character from windows
sed -i 's/\r//g' "$1"

# remove spaces after pipe character
sed -i 's/| /|/g' "$1"

# because of the next command it is important that the input file is of the shape file.fasta! if it's .fq or something strange stuff will happen
filename=$(basename "$1" .fasta)

# first awk removes al headers that have been seen before
# second awk changes all non nucleotide characters in the sequences
# third awk writes each different type of segment to a different segment-specific file (both header and segment)

awk '/^>/ {f=!d[$1]; d[$1]=1} f' $1 |  awk '/^>/ {print; next} {gsub(/[^ACGTURYSWKMBDHVN?-acgturyswkmbdhvn]/, "n"); print}'| awk -v filename="$filename" '/^>/{split($0, parts, "|"); segment = parts[length(parts)]; output_file =filename "_" segment ".fasta"; print $0 > output_file; next } {print > output_file}'
