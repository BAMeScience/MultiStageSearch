import pandas as pd
import numpy as np


peptide_shaker_report = snakemake.input[0]
strain_accessions = snakemake.input[1]

out_file = snakemake.output[0]

relevant_genomes = snakemake.params[0]


df = pd.read_csv(peptide_shaker_report, sep="\t", on_bad_lines="skip", usecols=["Protein(s)", "Confidence [%]"])
df.columns = ["ORF", "Confidence"]
df["weight"] = 1 / (df.ORF.str.count(",") + 1)
df["ORF"] = df.ORF.apply(lambda x: x.split(","))
df["psmid"] = df.index + 1
report_df = df.explode("ORF", ignore_index=True)
report_df["accession"] = report_df.ORF.apply(lambda x: x.split(".")[0])

genome_stats = pd.DataFrame({
    "accession": report_df["accession"],
    "counts": report_df.groupby("accession")["accession"].transform("count")
})


mean_conf_df = report_df.groupby("accession")["Confidence"].mean().reset_index()
mean_conf_df.rename(columns={"Confidence": "mean_confidence"}, inplace=True)
genome_stats = genome_stats.merge(mean_conf_df, on="accession")
genome_stats.drop_duplicates(inplace=True)
genome_stats.reset_index(drop=True, inplace=True)

accessions_df = pd.read_csv(strain_accessions, sep="\t")
accessions_df["genbank_accession"] = accessions_df.genbank_accession.apply(lambda x: x.split(".")[0])

merged_df = pd.merge(accessions_df, genome_stats, how="inner", left_on="genbank_accession", right_on="accession")
merged_df.replace(np.nan, "NaN", regex=True, inplace=True)
merged_df.sort_values("counts", ascending=False, inplace=True)
merged_df.reset_index(drop=True, inplace=True)
merged_df.drop(columns=["accession"])
merged_df = merged_df[:relevant_genomes]
merged_df.to_csv(out_file, index=False, sep="\t")