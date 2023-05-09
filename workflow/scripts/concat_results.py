import pandas as pd

firstDBSearch = snakemake.input[0]
finalDBSearch = snakemake.input[1]

concat_res = snakemake.output[0]


df1 = pd.read_csv(firstDBSearch, sep="\t", index_col=0)
df2 = pd.read_csv(finalDBSearch, sep="\t", index_col=0)
merged_df = pd.concat([df1, df2])
merged_df.reset_index(drop=True, inplace=True)

merged_df.to_csv(concat_res, sep="\t", header=True, index=True)