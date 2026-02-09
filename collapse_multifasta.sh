# !/bin/bash

# converting windows format to unix format by removing the the newline character from windows
sed -i 's/\r//g' "$1"

# remove spaces after pipe character
sed -i 's/| /|/g' "$1"

input=$1
output=$2

# first awk is to remove all those with same header
# sed removes the gene segment name of the header
# next awk concatanates the sequences with same header
# final awk is to replace all characters that aren't nucleotides or degenerated nucleotides

awk '/^>/ {f=!d[$1]; d[$1]=1} f' $input | sed -e 's/\(>.*\)|[^|]*$/\1/' | awk '/^>/ {if (!d[$1]++) print $0; next} !d[$1]' | awk '/^>/ {print; next} {gsub(/[^ACGTURYSWKMBDHVN?-acgturyswkmbdhvn]/, "n"); print}' > $outpu
