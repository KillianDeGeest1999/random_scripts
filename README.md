
script|goal|input|output|how to use
---|---|---|---|---
clean_fasta.sh|converts windows fasta to unix, removes spaces after pipe, removes doubles, removes all non nucleotides|"input_file".fasta|"input_file"_cleaned.fasta|bash clean_fasta.sh "input_file".fasta
cleaner.py|removes segment names from a fasta|"input_file".fasta|"input_file"_cleaned.fasta|python cleaner.py "input_file".fasta
collapse_multifasta.sh|windows to unix, removes spaces after pipe, removes doubles, concatenates segments in full genome, removes non nucleotides|"input_file".fasta|"output_file.fasta"|bash collapse_multifasta.sh "input_file".fasta "output_file".fasta
compare_fastas.py | compares sequences with matching IDs between two FASTA files and reports differences using pairwise alignment | file1.fasta + file2.fasta | sequence_differences.txt (alignment + length differences) | python compare_fastas.py
concatenate3.py|filters old samples if dates are included, orders the segments, concatenates the segments|"input_file".fasta|"output_file".fasta|run local
concatenator.py|concatenates segments in order|'input_file'.fasta|concatenated_'input_file'.fa|python concatenator.py 'input_file'.fasta
concatenator2.py|concatenates segments in order|'input_file'.fasta|concatenated_multifasta.fa|python concatenator.py 'input_file'.fasta
fasta_header_fixer.py|splits by "/" and replaces start with A, end with only first 3 characters, replaces all underscores with spaces|"input_file".fasta|"output_file".fasta|run local
remove_double.sh|converts windows fasta to unix, removes spaces after pipe, removes doubles|"input_file".fasta|"output_file".fasta|bash remove_double.sh "input_file".fasta "output".fasta
split_multifasta.sh|converts windows fasta to unix, removes spaces after pipe, removes doubles, removes all non nucleotides, splits multifasta per segment|"input_file".fasta|"input_file"_"segment".fasta|bash split_multifasta.sh "input_file".fasta
subsampler.py|removes highly similar sequences within same country-year-month groups (optionally keeps all Belgian sequences) based on pairwise alignment similarity threshold | 'input_file'.fasta similarity-score keep-belgium | 'output_file'.fasta| python subsampler.py input.fasta output.fasta 99 yes
subset_header.py | extracts selected fields from FASTA headers based on pipe-separated indices | 'input_file'.fasta header-indices | 'output_file'.fasta| python subset_header.py input.fasta output.fasta 0 2
trimmer.py | trims sequences to coding regions by keeping from first start codon (ATG) to last stop codon (TAA, TAG, TGA) | 'input_file'.fasta | 'output_file'.fasta | python trimmer.py 'input_file'.fasta 'output_file'.fasta
welke_mis_ik.py | identifies FASTA sequence headers that are missing from a reference list in an Excel file | 'input_file'.fasta + 'input_file'.xlsx | printed list of missing headers | python welke_mis_ik.py 'input_file'.fasta 'input_file'.xlsx
