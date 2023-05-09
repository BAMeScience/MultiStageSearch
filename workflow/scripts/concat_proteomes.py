import os
import shutil

concat_proteomes = snakemake.output[0]
res_dir = str(snakemake.params[0])
sample_name = str(snakemake.params[1])

cwd = str(os.getcwd())

full_in_path = f"{cwd}/{res_dir}/{sample_name}/sixpack"

included_extensions = ["fasta"]
file_names = [
    fn for fn in os.listdir(full_in_path) if any(fn.endswith(ext) for ext in included_extensions)
]

with open(concat_proteomes, "w") as out_file:
    for file in file_names:
        with open(f"{full_in_path}/{file}", "r") as f:
            shutil.copyfileobj(f, out_file)
