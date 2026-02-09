from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Align
from itertools import combinations
import sys
import datetime

# for the header we expect the shape to be number|sample_info|date|segment
def extract_country_and_year(header):
    parts = header.split('|') # here I split the header line on every pipe character and create a list with the different parts
    sample_info = parts[1].split('/') # here I split the part with the sample_ID for every /
    if len(sample_info)== 1: # there are a number of samples that are named wrongly and don't have / but, for example, underscores. These can't be split correctly and thus have a length of 1.
        country = parts[1] # As it is not possible to know how they f*cked it up I just save the entire sample id as country
    else:
        country = sample_info[2] # if normally split the country should be found on the third spot within the sample info.
    year = parts[2].split('-')[0]  # Extract year from the sample ID
    month = parts[2].split('-')[1]
    return country, year, month

def calculate_similarity(seq1, seq2): # this is a function that calculates the homology between two sequences and returns a score out of a 100
    aligner = Align.PairwiseAligner()
    alignments = aligner.align(seq1, seq2)
    best_score =alignments.score
    similarity = best_score / max(len(seq1), len(seq2)) * 100
    return similarity

def remove_similar_sequences(input_file, output_file, similarity_threshold, keep_Belgium):
    # ct stores current time
    ct = datetime.datetime.now()
    print("current time:-", ct) # info 
    records = list(SeqIO.parse(input_file, "fasta"))
    number_of_sequences=len(records)
    print(f' there are {number_of_sequences} sequences in the input file') # just as info
    unique_samples = set()
    number_removed=0

    # Extract country and year information for each sequence
    for record in records:
        country, year, month = extract_country_and_year(record.description)
        unique_samples.add((country, year, month))

    # Iterate through each unique combination of country and year
    indices_to_remove = set()
    print(keep_Belgium)
    if keep_Belgium == 'yes' : # if we want to keep all of the Belgian sequences we don't compare them. There should definitely be a better way to do this but it works
        for country, year, month in unique_samples:
            if country not in ['Belgium']:
                seq_indices = [i for i, record in enumerate(records) if extract_country_and_year(record.description) == (country, year, month)]
                for i, j in combinations(seq_indices, 2):
                    seq1 = records[i].seq
                    seq2 = records[j].seq
                    if i not in indices_to_remove and j not in indices_to_remove:
                        similarity = calculate_similarity(seq1, seq2)
                        if similarity >= float(similarity_threshold):
                            if len(seq1) >= len(seq2):
                                indices_to_remove.add(j)
                            else:
                                indices_to_remove.add(i)
                            number_removed+=1
                            if number_removed % 1000 == 0:
                                print(number_removed)
    elif keep_Belgium =='no': 
        for country, year, month in unique_samples:
            seq_indices = [i for i, record in enumerate(records) if extract_country_and_year(record.description) == (country, year, month)]
            for i, j in combinations(seq_indices, 2):
                seq1 = records[i].seq
                seq2 = records[j].seq
                if i not in indices_to_remove and j not in indices_to_remove:
                    similarity = calculate_similarity(seq1, seq2)
                    if similarity >= float(similarity_threshold):
                        if len(seq1) >= len(seq2):
                            indices_to_remove.add(j)
                        else:
                            indices_to_remove.add(i)
                        number_removed+=1
                        if number_removed % 1000 == 0:
                            print(number_removed)
    print(f'We have removed {number_removed} sequences that had a similarity above {similarity_threshold}, this is {(number_removed/number_of_sequences)*100} percent of the sequences.')

    # Write non-similar sequences to output file
    with open(output_file, "w") as output_handle:
        for i, record in enumerate(records):
            if i not in indices_to_remove:
                SeqIO.write(record, output_handle, "fasta")
    # ct stores current time
    ct = datetime.datetime.now()
    print("current time:-", ct)  
       
if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python script.py input.fasta output.fasta threshold_similarity keep_Belgium (keep_Belgium should be a yes or no)" )
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    similarity_threshold = sys.argv[3]
    keep_Belgium = sys.argv[4]
    remove_similar_sequences(input_file, output_file, similarity_threshold, keep_Belgium)
