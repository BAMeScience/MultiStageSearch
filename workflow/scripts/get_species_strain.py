import pandas as pd
from Bio import Entrez, SeqIO


in_file = snakemake.input[0]
out_file = snakemake.output[0]

relevant_genomes = snakemake.params[0]
APImail = snakemake.params[1]
APIkey = snakemake.params[2]


if APIkey and APImail:
    Entrez.email = APImail
    Entrez.api_key = APIkey


df = pd.read_csv(in_file, sep="\t", on_bad_lines="skip", usecols=["Protein(s)"])
df.columns = ["ORF"]
df["weight"] = 1 / (df.ORF.str.count(",") + 1)
df["ORF"] = df.ORF.apply(lambda x: x.split(","))
df["psmid"] = df.index + 1
report_df = df.explode("ORF", ignore_index=True)
report_df["accession"] = report_df.ORF.apply(lambda x: x.split(".")[0])

genome_mapping = pd.DataFrame()
genome_mapping["counts"] = report_df["accession"].value_counts()[:relevant_genomes]
genome_mapping = genome_mapping.reset_index()
genome_mapping = genome_mapping.rename(columns={"index" : "genbank_accession"})

new_df = pd.DataFrame(columns=["genbank_accession", "species", "strain", "isolate"])


for accession in genome_mapping["genbank_accession"]:
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")

    record = SeqIO.read(handle, "genbank")

    # Extract the strain name from the source feature
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

merged_df = pd.merge(genome_mapping, new_df, how="inner", left_on="genbank_accession", right_on="genbank_accession")

merged_df.to_csv(out_file, index=False, sep="\t")