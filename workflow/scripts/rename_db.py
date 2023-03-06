import os


res_dir = str(snakemake.params[0])
sample_name = str(snakemake.params[1]) 

cwd = str(os.getcwd())

full_in_path = f"{cwd}/resources/Database/{sample_name}"
full_out_path = f"{cwd}/{res_dir}/{sample_name}/Database"

included_extensions = ['fasta']
file_names = [fn for fn in os.listdir(full_in_path)
              if any(fn.endswith(ext) for ext in included_extensions)]

os.rename(full_in_path + "/" + file_names[0], full_out_path + "/" + sample_name + "_protein_concatenated_target_decoy.fasta")
