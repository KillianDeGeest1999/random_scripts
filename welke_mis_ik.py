import pandas as pd
from Bio import SeqIO
import sys

def main(fasta_file, excel_file):
    # Step 1: Parse headers from the multi-FASTA file
    fasta_headers = [record.id for record in SeqIO.parse(fasta_file, 'fasta')]

    # Step 2: Read the first column of the Excel file
    excel_data = pd.read_excel(excel_file, usecols=[0])  # Reads only the first column
    # Replace spaces with underscores in the Excel headers and strip '>' from the start
    excel_headers = excel_data.iloc[:, 0].astype(str).str.replace(' ', '_').str.lstrip('>').tolist()

    # Step 3: Compare headers and find those not in the Excel file
    missing_headers = [header for header in fasta_headers if header not in excel_headers]

    # Step 4: Print the missing headers
    print("Headers from the multi-FASTA file not found in the first column of the Excel file:")
    for header in missing_headers:
        print(header)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python check_headers.py <path_to_fasta_file> <path_to_excel_file>")
        sys.exit(1)

    fasta_file = sys.argv[1]
    excel_file = sys.argv[2]
    
    main(fasta_file, excel_file)
