import os

mgf_file = snakemake.input[0]
concat_report = snakemake.input[1]
ref = snakemake.input[2]
par_file = snakemake.input[3]
host_name = snakemake.params[0]
sample_name = snakemake.params[1]

cwd = os.getcwd()
res_dir = snakemake.params[2]

out_file = snakemake.output[0]

res_dir = f"{cwd}/{res_dir}"

ref_name = ref.split("/")[-1]
ref_name = ref_name.split(".")
ref_name = "".join(ref_name[:-1])

sci_host_name = f'"{host_name}"'

#print(sample_name, host_name)
print(ref_name)

with open(out_file, "w") as out_f:
    out_f.write("APImail: '' # provide API mail (optional but recommended)\n")
    out_f.write("APIkey: '' # provide API key (optional but recommended)\n")
    out_f.write(f"ExperimentName: '{sample_name}'\n")
    out_f.write(f"SampleName: '{sample_name}'\n")
    out_f.write(f"HostName: '{host_name}'\n")
    out_f.write(f"ReferenceDBName: '{ref_name}'\n")
    out_f.write(f"ScientificHostName: '{sci_host_name}'\n")
    out_f.write("FilterSpectra: False\n")
    out_f.write("AddHostandCrapToDB: True\n\n")
    out_f.write(f"SamplePath: '{mgf_file}'\n")
    out_f.write(f"ParametersFile: '{par_file}'\n")
    out_f.write("DataDir: 'resources/SampleData/'\n")  # TODO
    out_f.write("DatabaseDir: 'resources/Database/'\n")
    out_f.write("PeptideShaker: '/path/to/bin/PeptideShaker-2.2.16/PeptideShaker-2.2.16.jar'\n")
    out_f.write("SearchGUI: '/path/to/bin/SearchGUI-4.1.23/SearchGUI-4.1.23.jar'\n")
    out_f.write("ResourcesDir: 'resources/'\n")
    out_f.write("ResultsDir: 'results/'\n")
    out_f.write("TaxidMapping: 'results/taxidMapping'\n")
    out_f.write("searchengines: '-xtandem'\n")
    out_f.write("psmFDR: '5.0'\n")
    out_f.write("peptideFDR: '5.0'\n")
    out_f.write("proteinFDR: '5.0'\n")
    out_f.write("TaxaInPlot: 15\n")
    out_f.write("Alpha: [0.01,0.05,0.10,0.20,0.40,0.60]\n")
    out_f.write("Beta: [0.01,0.05,0.10,0.20,0.40,0.50,0.70]\n")
    out_f.write("prior: [0.10,0.30,0.50]")
