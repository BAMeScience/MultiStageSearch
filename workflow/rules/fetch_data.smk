rule fetchGenomes:
    input: 
        top_scoring_taxIDs = RESULT_DIR / "{sample}/taxids/top_scoring_taxids.tsv",
    output:
        concat_fasta = RESULT_DIR / "{sample}/Database/concat_taxid_genomes.fasta"
    log:
        stderr_log=RESULT_DIR / "logs/FetchData/{sample}/stderr.log",
        stdout_log=RESULT_DIR / "logs/FetchData/{sample}/stdout.log"
    # params:
    #     number_taxids = config["taxids"]["number_of_taxids"]
    conda:
        "../envs/fetch_data.yml"
    threads: 1
    script:
        "../scripts/fetch_data.py"