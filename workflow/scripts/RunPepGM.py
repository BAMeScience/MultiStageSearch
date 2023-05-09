import os
import subprocess

in_file = snakemake.input[0]
#in_file = "results/sample/concatRes/concat_PSM_Report.txt"
out_file = snakemake.output[1]

pepgm_path = snakemake.params[0]
res_dir = snakemake.params[1]

sample_name = in_file.split("/")[1]
print(sample_name)

with open(out_file, "w") as f:
    os.chdir(pepgm_path)
    process = subprocess.Popen(
        [
            "snakemake",
            "-c",
            "--use-conda",
            "-n",
            "--config",
            f"samplename={sample_name}",
            "hostname=host",
            "DBname=refSeqViral",
            f"ResultsDir={res_dir}"
        ],
        stdout=f,
        stderr=subprocess.STDOUT
    )


# os.chdir('/tmp')

# """
# samplename
# hostname
# DBname
# ResultsDir
# NumberofResults
# """