import subprocess
import time

concat_fasta = snakemake.input[0]
sixpack_temp = snakemake.output[1]

sample_path = snakemake.params[0] + "/sixpack/"
orfminsize = snakemake.params[1]
additional_sixpack_params = snakemake.params[2]


with open(concat_fasta, "r") as fasta:
    sequences = fasta.read()
    d = ">"
    sequences = [d+e for e in sequences.split(d)][1:]

with open(sixpack_temp, "w") as f:
    for i, sequence in enumerate(sequences):
        if additional_sixpack_params:
            process = subprocess.Popen(
                [
                    "sixpack",
                    "-sequence",
                    concat_fasta,
                    "-outfile",
                    f"{sample_path}{sequence[1:12]}.orfs",
                    "-outseq",
                    f"{sample_path}{sequence[1:12]}.fasta",
                    "-orfminsize",
                    f"{orfminsize}",
                    f"{additional_sixpack_params}"
                ],
                stdout=f,
                stderr=subprocess.STDOUT
            )
        else:
            process = subprocess.Popen(
                [
                    "sixpack",
                    "-sequence",
                    concat_fasta,
                    "-outfile",
                    f"{sample_path}{sequence[1:12]}.orfs",
                    "-outseq",
                    f"{sample_path}{sequence[1:12]}.fasta",
                    "-orfminsize",
                    f"{orfminsize}",
                ],
                stdout=f,
                stderr=subprocess.STDOUT
            )
        with open(concat_fasta, "w") as fasta:
            for seq in sequences[i:]:
                fasta.write(seq + "\n")
        time.sleep(1)  # because file system latency
