from Bio import Entrez
from Bio import SeqIO
import pandas as pd
import time


accessions_list = snakemake.input[0]
genbank = snakemake.input[1]
mapped_taxids = snakemake.input[2]

number_taxids = snakemake.params[0]
weight_diff = snakemake.params[1]

APImail = snakemake.params[2]
APIkey = snakemake.params[3]

concat_fasta = snakemake.output[0]

if APIkey and APImail:
    Entrez.email = APImail
    Entrez.api_key = APIkey


genbank_accessions = pd.read_csv(genbank, sep="\t", header=0, names=["accession", "name", "others"], index_col=False)
genbank_accessions = genbank_accessions["accession"].to_list()

accessions_to_query = []
with open(accessions_list, "r") as f:
    accessions = f.readlines()
    for accession in accessions:
        accession = accession.strip()
        if accession in genbank_accessions:
            accessions_to_query.append(accession)

print(len(accessions_to_query))
print(accessions_to_query)

sequences = Entrez.efetch(db="nucleotide", id=accessions_to_query, rettype="fasta", retmax=1000).read()
num_seq = sequences.count(">")
print(f"There were {num_seq} genome sequences found.")

#TODO code auslagern da doppelt
if num_seq == 0:
    print("Trying to use the species level genomes!")
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
    for taxid in taxids_to_query:
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
