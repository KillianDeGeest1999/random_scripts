import sys
import re

def clean_headers(input_file):
    output_file = input_file.rsplit('.', 1)[0] + '_cleaned.fasta'

    with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
        for line in fin:
            if line.startswith('>'):
                # Replace specified patterns in the header line
                line = re.sub(r'\|PB1\||\|PB2\||\|NA\||\|HA\||\|NP\||\|MP\||\|NS\|', '|', line)
            fout.write(line)

    print("Output written to:", output_file)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py input_file.fasta")
        sys.exit(1)
    input_file = sys.argv[1]
    clean_headers(input_file)
