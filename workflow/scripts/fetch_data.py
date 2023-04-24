import pandas as pd
from Bio import Entrez, SeqIO


mapped_taxids = snakemake.input[0]
out_file_df = snakemake.output[0]
out_file_fasta = snakemake.output[1]
APImail = snakemake.params[0]
APIkey = snakemake.params[1]
number_taxids = snakemake.params[2]
weight_diff = snakemake.params[3]
sequence_length_diff = snakemake.params[4]
max_number_accessions = snakemake.params[5]


#DEBUG
# mapped_taxids = "/home/jpipart/project/MultiStageSearch/results/PXD003013_Cowpox_BR/taxids/taxid_scores.tsv"
# out_file_df = "/home/jpipart/project/MultiStageSearch/test_strain_df.tsv"
# out_file_fasta = "/home/jpipart/project/MultiStageSearch/test_concat.fasta"
# relevant_genomes = 10
# APImail = "julian.pipart@fu-berlin.de"
# APIkey = "8773c6c2cff7ccdcdda9dc42f4b561516909"
# number_taxids = 5
# weight_diff = 2
# sequence_length_diff = 2
# max_number_accessions = 200


if APIkey and APImail:
    Entrez.email = APImail
    Entrez.api_key = APIkey




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





all_records = []
#all_genbank_accessions = []
new_df = pd.DataFrame(columns=["genbank_accession", "species", "strain", "isolate"])

for taxid in taxids_to_query:
    search_query = f"txid{taxid}[Organism]"
    search_handle = Entrez.esearch(db="nucleotide", term=search_query, retmax=max_number_accessions, sort="Sequence Length")
    taxon_id_list = Entrez.read(search_handle)["IdList"]
    
    for id in range(len(taxon_id_list)):
        print(f"{id+1}/max{len(taxon_id_list)}")
        handle = Entrez.efetch(db="nucleotide", id=taxon_id_list[id], rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        if (id > 0) and (len(all_records[id-1].seq) >= (len(record.seq) * sequence_length_diff)):
            print(len(all_records[id-1].seq), (len(record.seq) * sequence_length_diff))
            print("Abort query here since the difference in sequence lenght is too big!")
            break
        all_records.append(record)


# get df with strain names...
for record in all_records:
    #all_genbank_accessions.append(record.id)
    accession = record.id
    species = record.features[0].qualifiers["organism"][0]
    try:
        strain = record.features[0].qualifiers["strain"][0]
    except KeyError:
        strain = "NaN"
    try:
        isolate = record.features[0].qualifiers["isolate"][0]
    except KeyError:
        isolate = "NaN"
    new_df = new_df.append({"genbank_accession": accession, "species": species, "strain": strain, "isolate": isolate}, ignore_index=True)


# get unique strains...
unique_df = new_df.drop_duplicates(subset=["species", "strain", "isolate"], keep="first").reset_index(drop=True)
unique_strain_records = [record for record in all_records if (record.id in unique_df["genbank_accession"].to_list())]

# write files...
unique_df.to_csv(out_file_df, index=False, sep="\t")
SeqIO.write(unique_strain_records, out_file_fasta, "fasta")