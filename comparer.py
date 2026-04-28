from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

file1 = {rec.id: str(rec.seq) for rec in SeqIO.parse("20260213.fa", "fasta")}
file2 = {rec.id: str(rec.seq) for rec in SeqIO.parse("20260213_old.fa", "fasta")}

with open("sequence_differences.txt", "w") as out:
    for k in file1:
        seq1 = file1[k]
        seq2 = file2.get(k)

        if seq2 is None:
            out.write(f"\n{k} missing in second file\n")
            continue

        if seq1 != seq2:
            out.write(f"\n===== Difference in {k} =====\n")
            out.write(f"Length file1: {len(seq1)}\n")
            out.write(f"Length file2: {len(seq2)}\n\n")

            alignments = pairwise2.align.globalxx(seq1, seq2)
            out.write(format_alignment(*alignments[0]))
            out.write("\n")

print("Done. Results written to sequence_differences.txt")
