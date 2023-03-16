from Bio import Entrez
import pandas as pd
import sys

mapped_taxids = snakemake.input[0]
concat_fasta = snakemake.output[0]
number_taxids = snakemake.params[0]
stderr_log = snakemake.log[0]

with open(stderr_log, "w") as sys.stderr:

    taxid_df = pd.read_csv(mapped_taxids, sep="\t", header=0, index_col=0)
    taxids = taxid_df["taxid"].value_counts()[:number_taxids]
    taxids = taxids.index.to_list()

    genome_ids = []
    for taxid in taxids:
        q = f"refseq[filter] AND txid{taxid}[ORGN]"
        handle = Entrez.esearch(db="nucleotide", term=q, retmax=1000)
        genome_ids.append(int(handle.readlines()[3][4:-6]))
        handle.close()

    sequences = Entrez.efetch(db="nucleotide", id=genome_ids, rettype="fasta", retmax=1000).read()
    with open(concat_fasta, "w") as f:
        f.write(sequences)



# def score(df, taxids, output_path):
#     """
#     Score taxids according to their confidence and select the ones which are top scoring.
#     :param df: df, contains protein 2 taxid mapping and the respective weights
#     :param taxids: lst, mapped taxids
#     :param output_path: str, output path
#     :param subset: int, top <subset> scoring taxids
#     """
#     df['taxid'] = taxids
#     df_score = df.groupby('taxid')['weight'].sum().reset_index()
#     df_score = df_score.sort_values(by=['weight'], ascending=False)
#     threshold = df_score.weight.median()
#     top_scoring_df = df_score.loc[df_score["weight"] >= threshold]
#     top_scoring_df.to_csv(output_path[:-4] + '_weights.csv', index=False)
#     top_scoring_taxids = top_scoring_df.taxid.tolist()
#     save(output_path, top_scoring_taxids)