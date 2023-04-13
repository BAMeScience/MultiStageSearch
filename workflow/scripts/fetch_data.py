from Bio import Entrez
from Bio import SeqIO
import pandas as pd
import sys
import time

mapped_taxids = snakemake.input[0]
concat_fasta = snakemake.output[0]
number_taxids = snakemake.params[0]
weight_diff = snakemake.params[1]
stderr_log = snakemake.log[0]
stdout_log = snakemake.log[1]

# mapped_taxids = "/home/jpipart/project/MultiStageSearch/results/PXD025130_Sars_CoV_2/taxids/top_scoring_taxids.tsv"
# number_taxids = 5

# with open(stderr_log, "w") as sys.stderr:
#     with open(stdout_log, "w") as sys.stdout:

taxid_df = pd.read_csv(mapped_taxids, sep="\t", header=0, index_col=False)
#taxids = taxid_df["taxid"].value_counts()[:number_taxids]
taxids = taxid_df["taxid"].to_list()
relevant_taxids = taxids[:number_taxids]
taxids_to_query = []

for i in range(len(relevant_taxids)):
    if i == 0:
        taxids_to_query.append(relevant_taxids[i])
    else:
        previous_taxid = taxid_df.loc[taxid_df["taxid"] == relevant_taxids[i-1]]
        previous_taxid_score = previous_taxid["weight"].values[0]
        current_taxid = taxid_df.loc[taxid_df["taxid"] == relevant_taxids[i]]
        current_taxid_score = current_taxid["weight"].values[0]

        if previous_taxid_score <= (current_taxid_score * weight_diff):
            taxids_to_query.append(relevant_taxids[i])
        else:
            break

print("TaxIDs to query: ", taxids_to_query)


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
