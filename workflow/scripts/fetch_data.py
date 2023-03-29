from Bio import Entrez
from Bio import SeqIO
import pandas as pd
import sys
import time

mapped_taxids = snakemake.input[0]
concat_fasta = snakemake.output[0]
#number_taxids = snakemake.params[0]
stderr_log = snakemake.log[0]
stdout_log = snakemake.log[1]

# with open(stderr_log, "w") as sys.stderr:
#     with open(stdout_log, "w") as sys.stdout:

taxid_df = pd.read_csv(mapped_taxids, sep="\t", header=0, index_col=False)
#taxids = taxid_df["taxid"].value_counts()[:number_taxids]
taxids = taxid_df["taxid"].to_list()

genome_ids = []
for taxid in taxids:
    q = f"refseq[filter] AND txid{taxid}[ORGN]"
    handle = Entrez.esearch(db="nucleotide", term=q, retmax=1000)
    complete_handle = handle.readlines()
    if "<Count>0" in str(complete_handle[2]):
        print(f"No genome found for taxid {taxid}")
        continue
    genome_id = int(complete_handle[3][4:-6])
    genome_ids.append(genome_id)
    handle.close()

sequences = Entrez.efetch(db="nucleotide", id=genome_ids, rettype="fasta", retmax=1000).read()
num_seq = sequences.count(">")
print(f"There were {num_seq} genome sequences found.")
with open(concat_fasta, "w") as f:
    f.write(sequences)


# filter if there are accidental duplicates

records = list(SeqIO.parse(concat_fasta, "fasta"))
seen = set()
filtered_seqs = []

for record in records:  
    if record.name not in seen:
        seen.add(record.name)
        filtered_seqs.append(record)

SeqIO.write(filtered_seqs, concat_fasta, "fasta")

# TODO: Need a better solution for this
with open(concat_fasta, "r") as fasta:
    sequences = fasta.read()
    d = ">"
    s = [d+e for e in sequences.split(d)]

with open(concat_fasta, "w") as fasta:
    for seq in s[1:]:
        # print(seq[:20])
        fasta.write(seq)
time.sleep(1)
