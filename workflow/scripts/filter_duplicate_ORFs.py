from Bio import SeqIO

in_file = snakemake.input[0]
out_file = snakemake.output[0]

sequences = {}
counter = 0
for record in SeqIO.parse(in_file, "fasta"):
    counter += 1
    if record.seq not in sequences:
        sequences[record.seq] = record.id
    else:
        sequences[record.seq] += f",{record.id}"

print(f"Number of sequences: {counter}")
print(f"Number of unique sequences: {len(sequences)}")


# write to fasta file
with open(out_file, "w") as f:
    for sequence in sequences:
        seq = str(sequence)
        f.write(">" + str(sequences[sequence]) + "\n")
        while seq:
            f.write(seq[:60] + "\n")
            seq = seq[60:]
        f.write("\n")