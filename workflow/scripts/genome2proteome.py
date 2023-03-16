import subprocess
import time

# concat_fasta = "/home/jpipart/project/MultiStageSearch/Database/concat_taxid_genomes.fasta"
# sixpack_temp = "sixpack_temp.txt"

# base_path = "results/PXD003013_Cowpox_BR/sixpack"
# sample = ""

concat_fasta = snakemake.input[0]
sixpack_temp = snakemake.output[1]

sample_path = snakemake.params[0] + "/sixpack/"


with open(concat_fasta, "r") as fasta:
    sequences = fasta.read()
    sequences = sequences.split("\n\n")


with open(sixpack_temp, "w") as f:
    for i, sequence in enumerate(sequences[:-1]):
        process = subprocess.Popen(
            [
                "sixpack",
                "-sequence",
                concat_fasta,
                "-outfile",
                f"{sample_path}{sequence[1:12]}.orfs",
                "-outseq",
                f"{sample_path}{sequence[1:12]}.fasta",
            ],
            stdout=f,
            stderr=subprocess.STDOUT
        )
        with open(concat_fasta, "w") as fasta:
            for seq in sequences[i:]:
                # print(seq[:20])
                fasta.write(seq + "\n\n")
        time.sleep(1)
