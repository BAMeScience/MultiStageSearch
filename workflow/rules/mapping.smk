rule mapTaxIDs:
    input:
        PSM_Report = RESULT_DIR / "{sample}/FirstSearch/ref_Default_PSM_Report.txt",
        taxIDs = RESULT_DIR / "taxidMapping/protacc2taxids_virus.txt"
    output:
        all_mapped_taxIDs = RESULT_DIR / "{sample}/taxids/mapped_taxids.tsv",
        top_scoring_taxIDs = RESULT_DIR / "{sample}/taxids/top_scoring_taxids.tsv",
    log:
        stderr_log = RESULT_DIR / "logs/taxids/getTargets/{sample}/stderr.log",
        stdout_log = RESULT_DIR / "logs/taxids/getTargets/{sample}/stdout.log"
    params:
        number_taxids = config["mapping"]["number_of_taxids"]
    conda: 
        "../envs/mapping.yml"
    script: 
        "../scripts/taxIDMapping.py"