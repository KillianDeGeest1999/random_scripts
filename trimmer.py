from Bio import SeqIO
import sys


def find_final_stop_codon(sequence):
    stop_codons = ["TAA", "TAG", "TGA"]
    final_stop = -1
    for codon in stop_codons:
        stop_pos = sequence.upper().rfind(codon)
        if stop_pos != -1 and (final_stop == -1 or stop_pos > final_stop):
            final_stop = stop_pos
    return final_stop

def trim_sequences(input_file, output_file):
    with open(output_file, "w") as output_handle:
        for record in SeqIO.parse(input_file, "fasta"):
            sequence = record.seq
            sequence_upper = sequence.upper()
            
            # Find start codon
            start_codon = sequence_upper.find("ATG")
            
            # Find final stop codon
            final_stop = find_final_stop_codon(sequence_upper)
            
            # Trim the sequence if both start and stop codons are found
            if start_codon != -1 and final_stop != -1:
                trimmed_sequence = sequence[start_codon:final_stop+3]
                record.seq = trimmed_sequence
            SeqIO.write(record, output_handle, "fasta")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py input.fasta output.fasta")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    trim_sequences(input_file, output_file)
