rule CreatePepGMConfig:
    input:
        mgf = MGF_FILE,
        concatResults = RESULT_DIR / "{sample}/concatRes/concat_PSM_Report.txt",
        ref = config["db_search"]["ref"],
        par = PAR_FILE,
    output:
        saple_config = RESULT_DIR / "{sample}/PepGM/config/config.yaml"
    params:
        host_name = HOST_NAME,
        sample_name = "{sample}",
        res_dir = str(RESULT_DIR)
    script:
        "../scripts/createPepGMConfig.py"




from snakemake.utils import min_version
min_version("6.0")

# module PepGM:
#     snakefile:
#         # here, plain paths, URLs and the special markers for code hosting providers (see below) are possible.
#         "/home/jpipart/PepGM/workflow/Snakefile"

# use rule CreateFactorGraph from PepGM as PepGM_CreateFactorGraph



module PepGM:
    snakefile:
        #github("BAMeScience/PepGM", path="workflow/Snakefile", branch="master")#, tag="v2.0.1")
         "/home/jpipart/PepGM/workflow/Snakefile"
    config: RESULT_DIR / "{sample}/PepGM/config/config.yaml"
    prefix: "PepGM"


use rule CreateFactorGraph from PepGM as PepGM_CreateFactorGraph with:
    input:
        RESULT_DIR / "{sample}/refSeqViral_Default_PSM_Report.txt",
        RESULT_DIR / "{sample}/taxids/mapped_taxids.tsv"
    output:
        RESULT_DIR / "{sample}/PepGM/PepGM_graph.graphml"


# rule CreateFactorGraph:
#     input: ResultsDir + SampleName + '/{DBname}_Default_PSM_Report.txt',
#            ResultsDir + SampleName + '/{DBname}_mapped_taxids.txt'
#     output: ResultsDir + SampleName + '/{DBname}_PepGM_graph.graphml'


# rule RunPepGM:
#     input:
#         concatResults = RESULT_DIR / "{sample}/concatRes/concat_PSM_Report.txt"
#     output:
#         touch(RESULT_DIR / "{sample}/PepGM/PepGM_run.done"),
#         PepGM_run = RESULT_DIR / "{sample}/PepGM/PepGM_run.txt"
#     params:
#         PepGM_path = config["pepgm_path"],
#         res_dir = RESULT_DIR / "{sample}/PepGM/"
#     script:
#         "../scripts/RunPepGM.py"
