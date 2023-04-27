import pandas as pd
import numpy as np


peptide_shaker_report = snakemake.input[0]
strain_accessions = snakemake.input[1]

out_file = snakemake.output[0]

relevant_genomes = snakemake.params[0]



#DEBUG
# peptide_shaker_report = "/home/jpipart/project/MultiStageSearch/results/PXD002936_avian_bronchitis/FinalSearch/proteomes_Default_PSM_Report.txt"
# strain_accessions = "/home/jpipart/project/MultiStageSearch/results/PXD002936_avian_bronchitis/FetchData/strain_accessions.tsv"
# out_file = "/home/jpipart/project/MultiStageSearch/test_strain_names.tsv"
# relevant_genomes = 20


df = pd.read_csv(peptide_shaker_report, sep="\t", on_bad_lines="skip", usecols=["Protein(s)"])
df.columns = ["ORF"]
df["weight"] = 1 / (df.ORF.str.count(",") + 1)
df["ORF"] = df.ORF.apply(lambda x: x.split(","))
df["psmid"] = df.index + 1
report_df = df.explode("ORF", ignore_index=True)
report_df["accession"] = report_df.ORF.apply(lambda x: x.split(".")[0])

genome_counts = pd.DataFrame()
genome_counts["counts"] = report_df["accession"].value_counts()[:relevant_genomes]
genome_counts = genome_counts.reset_index()
genome_counts = genome_counts.rename(columns={"index" : "genbank_accession"})

accessions_df = pd.read_csv(strain_accessions, sep="\t")#, on_bad_lines="skip")

accessions_df["genbank_accession"] = accessions_df.genbank_accession.apply(lambda x: x.split(".")[0])

merged_df = pd.merge(accessions_df, genome_counts, how="inner", left_on="genbank_accession", right_on="genbank_accession")
merged_df.replace(np.nan, "NaN", regex=True, inplace=True)
merged_df.sort_values("counts", ascending=False, inplace=True)
merged_df.reset_index(drop=True, inplace=True)
merged_df.to_csv(out_file, index=False, sep="\t")


# df = pd.read_csv(in_file, sep="\t", on_bad_lines="skip", usecols=["Protein(s)"])
# df.columns = ["ORF"]
# df["weight"] = 1 / (df.ORF.str.count(",") + 1)
# df["ORF"] = df.ORF.apply(lambda x: x.split(","))
# df["psmid"] = df.index + 1
# report_df = df.explode("ORF", ignore_index=True)
# report_df["accession"] = report_df.ORF.apply(lambda x: x.split(".")[0])

# genome_mapping = pd.DataFrame()
# genome_mapping["counts"] = report_df["accession"].value_counts()[:relevant_genomes]
# genome_mapping = genome_mapping.reset_index()
# genome_mapping = genome_mapping.rename(columns={"index" : "genbank_accession"})

# new_df = pd.DataFrame(columns=["genbank_accession", "species", "strain", "isolate"])


# for accession in genome_mapping["genbank_accession"]:
#     handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")

#     record = SeqIO.read(handle, "genbank")

#     # Extract the strain name from the source feature
#     species = record.features[0].qualifiers["organism"][0]
#     try:
#         strain = record.features[0].qualifiers["strain"][0]
#     except KeyError:
#         strain = "NaN"
#     try:
#         isolate = record.features[0].qualifiers["isolate"][0]
#     except KeyError:
#         isolate = "NaN"
#     new_df = new_df.append({"genbank_accession": accession, "species": species, "strain": strain, "isolate": isolate}, ignore_index=True)

# merged_df = pd.merge(genome_mapping, new_df, how="inner", left_on="genbank_accession", right_on="genbank_accession")

# merged_df.to_csv(out_file, index=False, sep="\t")