import pandas as pd

raw_report = snakemake.input[0]
taxids = snakemake.input[1]
out_file = snakemake.output[0]

raw_report_df = pd.read_csv(raw_report, sep = "\t", on_bad_lines="skip", usecols=["Protein(s)"])
id_df = pd.read_csv(taxids, header=None, names=["protein", "taxid"], sep="\t")

raw_report_df.columns = ["accession"]
raw_report_df["weight"] = 1 / (raw_report_df.accession.str.count(",") + 1)
raw_report_df["accession"] = raw_report_df.accession.apply(lambda x: x.split(","))
raw_report_df["psmid"] = raw_report_df.index + 1
report_df = raw_report_df.explode("accession", ignore_index=True)

merged_df = pd.merge(report_df, id_df, how="inner", left_on= "accession", right_on= "protein")
merged_df = merged_df.drop("protein", axis=1)

merged_df.to_csv(out_file, sep="\t", header=["accession", "weight", "psmid", "taxid"])

