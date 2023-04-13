rule fetchGenomes:
    input: 
        taxID_scores = RESULT_DIR / "{sample}/taxids/taxid_scores.tsv",
    output:
        concat_fasta = RESULT_DIR / "{sample}/Database/concat_taxid_genomes.fasta"
    log:
        stderr_log=RESULT_DIR / "logs/FetchData/{sample}/stderr.log",
        stdout_log=RESULT_DIR / "logs/FetchData/{sample}/stdout.log"
    params:
        number_taxids = config["mapping"]["number_of_taxids"],
        weight_diff = config["mapping"]["max_weight_differences"]
    conda:
        "../envs/fetch_data.yml"
    threads: 1
    script:
        "../scripts/fetch_data.py"