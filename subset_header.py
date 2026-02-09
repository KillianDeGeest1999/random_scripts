# "Usage: python script.py <input_file> <output_file> <index1> <index2> ... <indexN>"

import sys

def modify_fasta_headers(input_file, output_file, *indices):
    indices = [int(index) for index in indices]
    
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                header_parts = line.strip().split('|')
                new_header_parts = [header_parts[index] for index in indices if index < len(header_parts)]
                new_header = '|'.join(new_header_parts)
                if new_header_parts[0] == header_parts[0]:
                    outfile.write(new_header + '\n')
                else:
                    outfile.write('>' + new_header + '\n')
            else:
                outfile.write(line)

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python script.py <input_file> <output_file> <index1> <index2> ... <indexN>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    indices = sys.argv[3:]
    
    modify_fasta_headers(input_file, output_file, *indices)
