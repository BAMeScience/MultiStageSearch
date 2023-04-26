import sys
import pandas as pd
from ete3 import NCBITaxa
from Bio import Entrez, SeqIO


def fetchGenomes(all_records, taxon_id_list, sequence_length_diff):
    for id in range(len(taxon_id_list)):
            print(f"{id+1}/max{len(taxon_id_list)}")
            handle = Entrez.efetch(db="nucleotide", id=taxon_id_list[id], rettype="gb", retmode="text")
            record = SeqIO.read(handle, "genbank")
            if (id > 0) and (len(all_records[id-1].seq) >= (len(record.seq) * sequence_length_diff)):
                print(len(all_records[id-1].seq), (len(record.seq) * sequence_length_diff))
                print("Abort query here since the difference in sequence lengrg is too big!")
                break
            all_records.append(record)
    return all_records


def getTaxIdsToQuery(mapped_taxids, number_taxids, weight_diff):

    taxid_df = pd.read_csv(mapped_taxids, sep="\t", header=0, index_col=False)
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
    return taxids_to_query


def fetchData(taxids_to_query, max_number_accessions, sequence_length_diff):
    all_records = []
    for taxid in taxids_to_query:
        search_query = f"txid{taxid}[Organism]"
        search_handle = Entrez.esearch(db="nucleotide", term=search_query, retmax=max_number_accessions, sort="Sequence Length")
        taxon_id_list = Entrez.read(search_handle)["IdList"]
        
        all_records = fetchGenomes(all_records, taxon_id_list, sequence_length_diff)
    return all_records


def fetchDataNCBI(taxids_to_query, max_number_accessions, sequence_length_diff):
    ncbi = NCBITaxa()
    all_records = []
    
    for taxid in taxids_to_query:
        descendants = ncbi.get_descendant_taxa(taxid)
        descendant_names = ncbi.translate_to_names(descendants)
        
        for name in descendant_names:
            search_query = f"{name}"
            search_handle = Entrez.esearch(db="nucleotide", term=search_query, retmax=max_number_accessions, sort="Sequence Length")
            taxon_id_list = Entrez.read(search_handle)["IdList"]
            
            all_records = fetchGenomes(all_records, taxon_id_list, sequence_length_diff)
    return all_records


def getStrainNames(all_records):
    new_df = pd.DataFrame(columns=["genbank_accession", "species", "strain", "isolate"])
    # get df with strain names...
    for record in all_records:
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
    unique_df = new_df.drop_duplicates(subset=["genbank_accession", "species", "strain", "isolate"], keep="first").reset_index(drop=True)
    unique_strain_records = [record for record in all_records if (record.id in unique_df["genbank_accession"].to_list())]

    return unique_df, unique_strain_records


def writeFiles(unique_df, unique_strain_records, out_file_df, out_file_fasta):
    # filter problematic records
    record_to_filter = []
    for record in range(len(unique_strain_records)):
        try:
            unique_strain_records[record].seq[0]
        except:
            record_to_filter.append(record)

    record_to_filter.sort(reverse=True)
    for record in record_to_filter:
        print(f"removing problematic sequence record {unique_strain_records[record]}")
        unique_strain_records.pop(record)
    # write files...
    unique_df.to_csv(out_file_df, index=False, sep="\t")
    SeqIO.write(unique_strain_records, out_file_fasta, "fasta")


def main():
    stderr_log = snakemake.log[0]
    stdout_log = snakemake.log[1]
    with open(stderr_log, "w") as sys.stderr:
        with open(stdout_log, "w") as sys.stdout:
            mapped_taxids = snakemake.input[0]
            out_file_df = snakemake.output[0]
            out_file_fasta = snakemake.output[1]
            APImail = snakemake.params[0]
            APIkey = snakemake.params[1]
            number_taxids = snakemake.params[2]
            weight_diff = snakemake.params[3]
            sequence_length_diff = snakemake.params[4]
            max_number_accessions = snakemake.params[5]
            use_NCBI_Taxa = snakemake.params[6]

            if APIkey and APImail:
                Entrez.email = APImail
                Entrez.api_key = APIkey

            taxids_to_query = getTaxIdsToQuery(mapped_taxids, number_taxids, weight_diff)

            if use_NCBI_Taxa:
                all_records = fetchDataNCBI(taxids_to_query, max_number_accessions, sequence_length_diff)
            else:
                all_records = fetchData(taxids_to_query, max_number_accessions, sequence_length_diff)

            unique_df, unique_strain_records = getStrainNames(all_records)
            writeFiles(unique_df, unique_strain_records, out_file_df, out_file_fasta)


if __name__ == "__main__":
    main()