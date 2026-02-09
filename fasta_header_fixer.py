
def modify_header(header):
    parts = header.split('/')
    
    # Replace the first part with 'A' if it hasn't happened yet
    if not parts[0].startswith('A'):
        parts[0] = 'A'
    # Trim the last part after the first three numbers
    parts[-1] = parts[-1][:3]
    # Combine the parts back into a string
    modified_header = '/'.join(parts)
    return modified_header
def replace_underscore_with_space(header):
    return header.replace('_', ' ')
def process_fasta_headers(input_file, output_file):
    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        headers_written = 0
        for line in f_in:
            if line.startswith('>'):
                # Split the header line and extract the header without '>'
                header = line.strip()[1:]
                # Modify the header
                modified_header = modify_header(header)
                # Replace underscores with spaces
                modified_header_with_space = replace_underscore_with_space(modified_header)
                # Write the modified header to the output file
                f_out.write(modified_header + ' "' + modified_header_with_space + '" ')
                # Increment the counter for headers written
                headers_written += 1
                # If 50 headers have been written, start a new line
                if headers_written % 50 == 0:
                    f_out.write('\n')
        # Add a newline character at the end in case the total number of headers is not a multiple of 50
        if headers_written % 50 != 0:
            f_out.write('\n')
# Replace 'input.fasta' with the path to your input FASTA file
# Replace 'output.fasta' with the desired path for the output FASTA file
process_fasta_headers('1_PB2_pilot_clean.fasta', 'pilotdata_search_100.txt')

# Replace 'input.fasta' with the path to your input FASTA file
# Replace 'output.fasta' with the desired path for the output FASTA file
process_fasta_headers('1_PB2_pilot_clean.fasta', 'output.fasta')
