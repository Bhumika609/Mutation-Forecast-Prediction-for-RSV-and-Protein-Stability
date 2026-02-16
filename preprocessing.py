from Bio import SeqIO
from Bio.Seq import Seq

input_file = "sequence.fasta"
output_file = "cleaned_sequences.fasta"

min_len = 14000
max_len = 16000
max_n_percent = 0.10

records = list(SeqIO.parse(input_file, "fasta"))

print("Total sequences:", len(records))

filtered = []

for record in records:
    seq = str(record.seq).upper()
    length = len(seq)

    n_percent = seq.count("N") / length

    if (min_len <= length <= max_len and
        n_percent <= max_n_percent):

        record.seq = Seq(seq)
        filtered.append(record)

print("After quality filtering:", len(filtered))

# Remove exact duplicates
unique_sequences = {}
for record in filtered:
    seq = str(record.seq)
    if seq not in unique_sequences:
        unique_sequences[seq] = record

final_records = list(unique_sequences.values())

print("After removing exact duplicates:", len(final_records))

SeqIO.write(final_records, output_file, "fasta")

print("Cleaned file saved as:", output_file)