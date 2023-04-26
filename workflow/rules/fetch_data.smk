# rule getAccessionsToQuery:
#     input: 
#         taxID_scores = RESULT_DIR / "{sample}/taxids/taxid_scores.tsv",
#     output:
#         accessions_list = RESULT_DIR / "{sample}/FetchData/accessions_list.txt"
#     log:
#         stderr_log=RESULT_DIR / "logs/FetchData/{sample}/getAccessionsToQuery/stderr.log",
#         stdout_log=RESULT_DIR / "logs/FetchData/{sample}/getAccessionsToQuery/stdout.log"
#     params:
#         number_taxids = config["mapping"]["number_of_taxids"],
#         weight_diff = config["mapping"]["max_weight_differences"],
#         max_number_accessions = config["fetchData"]["max_number_accessions"],
#         APIMail = config["Entrez"]["APIMail"],
#         APIKey = config["Entrez"]["APIKey"],
#     conda:
#         "../envs/fetch_data.yml"
#     threads: 1
#     script:
#         "../scripts/get_accessions_to_query.py"

rule fetchStrainGenomes:
    input:
        taxID_scores = RESULT_DIR / "{sample}/taxids/taxid_scores.tsv",
    output:
        out_file_df = RESULT_DIR / "{sample}/FetchData/strain_accessions.tsv",
        concat_fasta = RESULT_DIR / "{sample}/FetchData/concat_strain_genomes.fasta"
    log:
        stderr_log=RESULT_DIR / "logs/FetchData/{sample}/fetchStrainGenomes/stderr.log",
        stdout_log=RESULT_DIR / "logs/FetchData/{sample}/fetchStrainGenomes/stdout.log"
    params:
        APIMail = config["Entrez"]["APIMail"],
        APIKey = config["Entrez"]["APIKey"],
        number_taxids = config["mapping"]["number_of_taxids"],
        weight_diff = config["mapping"]["max_weight_differences"],
        sequence_length_diff = config["fetchData"]["sequence_length_diff"],
        max_number_accessions = config["fetchData"]["max_number_accessions"],
        NCBI = config["fetchData"]["use_NCBI_Taxa"]
    conda:
        "../envs/fetch_data.yml"
    threads: 1
    script:
        "../scripts/fetch_data.py"


rule get_strain_names:
    input:
        peptide_shaker_report = RESULT_DIR / "{sample}/FinalSearch/proteomes_Default_PSM_Report.txt",
    output:
        accessions_mapping = RESULT_DIR / "{sample}/taxids/strain_name_mappings.tsv"
    log:
        stderr_log=RESULT_DIR / "logs/FetchData/{sample}/get_strain_names/stderr.log",
        stdout_log=RESULT_DIR / "logs/FetchData/{sample}/get_strain_names/stdout.log"
    params:
        number_taxids = config["mapping"]["number_of_strains"],
        APIMail = config["Entrez"]["APIMail"],
        APIKey = config["Entrez"]["APIKey"],
    conda:
        "../envs/fetch_data.yml"
    threads: 1
    script:
        "../scripts/get_species_strain.py"
