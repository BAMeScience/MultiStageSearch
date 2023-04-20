rule getAccessionsToQuery:
    input: 
        taxID_scores = RESULT_DIR / "{sample}/taxids/taxid_scores.tsv",
    output:
        accessions_list = RESULT_DIR / "{sample}/FetchData/accessions_list.txt"
    log:
        stderr_log=RESULT_DIR / "logs/FetchData/{sample}/getAccessionsToQuery/stderr.log",
        stdout_log=RESULT_DIR / "logs/FetchData/{sample}/getAccessionsToQuery/stdout.log"
    params:
        number_taxids = config["mapping"]["number_of_taxids"],
        weight_diff = config["mapping"]["max_weight_differences"],
        max_number_accessions = config["fetchData"]["max_number_accessions"],
        APIMail = config["Entrez"]["APIMail"],
        APIKey = config["Entrez"]["APIKey"],
    conda:
        "../envs/fetch_data.yml"
    threads: 1
    script:
        "../scripts/get_accessions_to_query.py"

rule fetchStrainGenomes:
    input:
        accessions_list = RESULT_DIR / "{sample}/FetchData/accessions_list.txt",
        genbank = config["genbank_accessions"],
        taxID_scores = RESULT_DIR / "{sample}/taxids/taxid_scores.tsv",
    output:
        concat_fasta = RESULT_DIR / "{sample}/FetchData/concat_taxid_genomes.fasta"
    log:
        stderr_log=RESULT_DIR / "logs/FetchData/{sample}/fetchStrainGenomes/stderr.log",
        stdout_log=RESULT_DIR / "logs/FetchData/{sample}/fetchStrainGenomes/stdout.log"
    params:
        number_taxids = config["mapping"]["number_of_taxids"],
        weight_diff = config["mapping"]["max_weight_differences"],
        APIMail = config["Entrez"]["APIMail"],
        APIKey = config["Entrez"]["APIKey"],
    conda:
        "../envs/fetch_data.yml"
    threads: 1
    script:
        "../scripts/fetch_strain_genomes.py"