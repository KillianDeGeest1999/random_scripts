#!/bin/bash

#converting windows format to unix format by removing the the newline character from windows
sed -i 's/\r//g' "$1"

# remove spaces after pipe character
sed -i 's/| /|/g' "$1"

input=$1
output=$2

# first awk is to remove all those with same header


awk '/^>/ {f=!d[$1]; d[$1]=1} f' $input > $output
