from Bio import Entrez
import pandas as pd

mapped_taxids = snakemake.input[0]
all_genbank_accessions = snakemake.output[0]
number_taxids = snakemake.params[0]
weight_diff = snakemake.params[1]
max_number_accessions = snakemake.params[2]
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

print("TaxIDs used: ", taxids_to_query)

accessions = []
for taxon_id in taxids_to_query:
# Use Entrez to search for the virus strains in Taxonomy database
    search_query = f"txid{taxon_id}[Organism]"
    search_handle = Entrez.esearch(db="nucleotide", term=search_query, retmax=max_number_accessions)

    # Get the list of matching taxon IDs
    taxon_id_list = Entrez.read(search_handle)["IdList"]
    #print(taxon_id_list)
    num_ids = len(taxon_id_list)

    # Get the full record for each taxon ID
    for taxon_id in range(num_ids):
        record_handle = Entrez.efetch(db="nucleotide", id=taxon_id_list[taxon_id], retmode="xml")
        record = Entrez.read(record_handle)[0]
        accession = record["GBSeq_primary-accession"]
        print(accession, f"{taxon_id + 1}/{num_ids}")
        accessions.append(accession)

with open(all_genbank_accessions, "w") as f:
    for accession in accessions:
        f.write(accession + "\n")