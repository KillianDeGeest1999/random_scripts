
import sys
def read_multifasta(file_name):
    """
    Reads a multifasta file and returns a dictionary containing headers as keys and corresponding sequences as values.
    """
    headers = {}
    current_header = None
    with open(file_name, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                current_header = line
                if current_header not in headers:
                    headers[current_header] = []
            else:
                if current_header in headers:
                    headers[current_header].append(line)
    return headers
def sort_multifasta(headers):
    """
    Sorts the multifasta sequences alphabetically based on the first part of the header and gene name.
    """
    # Sort headers alphabetically based on the first part of the header
    sorted_headers = sorted(headers.keys(), key=lambda x: x.split('|')[0])
    sorted_multifasta = {}
    for header in sorted_headers:
        sequence = headers[header]
        # Extract the prefix of the header (everything before the last '|') and the gene name
        header_prefix = '|'.join(header.split('|')[:-1])
        gene = header.split('|')[-1].split('/')[-1]
        # Create a nested dictionary with gene names as keys and header prefixes as subkeys
        if gene not in sorted_multifasta:
            sorted_multifasta[gene] = {}
        sorted_multifasta[gene][header_prefix] = sequence
    return sorted_multifasta
def concatenate_sequences(sorted_multifasta):
    """
    Concatenates sequences based on the order of genes and headers.
    """
    # Define the order of genes
    order = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]
    concatenated_multifasta = {}
    for gene in order:
        # Iterate over each gene in the specified order
        if gene in sorted_multifasta:
            for header_prefix, sequence in sorted_multifasta[gene].items():
                # Create a dictionary to store concatenated sequences for each header prefix
                if header_prefix not in concatenated_multifasta:
                    concatenated_multifasta[header_prefix] = []
                concatenated_multifasta[header_prefix].extend(sequence)
    return concatenated_multifasta
def write_multifasta(concatenated_multifasta, output_file):
    """
    Writes concatenated sequences to a new multifasta file.
    """
    with open(output_file, 'w') as file:
        for header_prefix, sequences in concatenated_multifasta.items():
            # Write each header prefix and its concatenated sequence to the output file
            file.write(f"{header_prefix}\n")
            file.write(''.join(sequences) + '\n')
def main():
    # Check if the correct number of command-line arguments are provided
    if len(sys.argv) != 2:
        print("Usage: python script.py input_file")
        sys.exit(1)
    input_file = sys.argv[1]  # Get the input file name from command line
    output_file = "concatenated_multifasta.fa"  # Define the output file name
    # Read multifasta file
    multifasta_headers = read_multifasta(input_file)
    # Sort multifasta sequences
    sorted_multifasta = sort_multifasta(multifasta_headers)
    # Concatenate sequences
    concatenated_multifasta = concatenate_sequences(sorted_multifasta)
    # Write concatenated sequences to a new file
    write_multifasta(concatenated_multifasta, output_file)
    print(f"Concatenated sequences saved to {output_file}")
if __name__ == "__main__":
    main()
